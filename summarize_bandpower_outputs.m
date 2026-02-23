function out = summarize_bandpower_outputs(input_folder, config_path)
%SUMMARIZE_BANDPOWER_OUTPUTS Generate analysis-friendly batch summary tables.
%
% Usage:
%   summarize_bandpower_outputs(input_folder)
%   summarize_bandpower_outputs(input_folder, config_path)
%
% This reads per-subject outputs under:
%   <input_folder>/bandpower_outputs/<subject_id>/...
%
% and writes merged summary CSVs into:
%   <input_folder>/bandpower_outputs/summary/
%
% Output struct 'out' contains written file paths.

if nargin < 2
    config_path = 'config.json';
end

cfg = pipeline.load_config(config_path);
fp_out = pipeline.get_output_root(input_folder, cfg);
fp_sum = pipeline.get_summary_dir(input_folder, cfg);

% Ensure summary-stage logs are written to summary/pipeline.log
% (Users often run summarize_bandpower_outputs separately after batch runs.)
didStartDiary = false;
try
    if isfield(cfg,'save_log') && logical(cfg.save_log)
        % Timestamped log file (never overwrite previous logs)
        ts = datestr(now, 'yyyymmdd_HHMMSS');
        logPathRun = fullfile(fp_sum, sprintf('pipeline_%s.log', ts));
        logPathLatest = fullfile(fp_sum, 'pipeline_latest.log');

        diaryOn = false;
        try
            diaryOn = strcmpi(get(0,'Diary'), 'on');
        catch
        end

        if ~diaryOn
            diary(logPathRun);
            diary on;
            didStartDiary = true;
        else
            % If diary is already on, do not switch files silently.
            try
                cur = string(get(0,'DiaryFile'));
                fprintf(2, '[NOTE] Diary is already ON: %s\n', cur);
                fprintf(2, '       This summarize run will NOT redirect diary to: %s\n', logPathRun);
            catch
            end
        end

        fprintf('\n===== summarize_bandpower_outputs started: %s =====\n', datestr(now,31));
        fprintf('Input folder: %s\n', input_folder);
        fprintf('Summary dir : %s\n', fp_sum);
        fprintf('Log file    : %s\n\n', logPathRun);
    end
catch
end

% List subject folders
D = dir(fp_out);
D = D([D.isdir]);
subs = setdiff({D.name}, {'.','..','summary'});
subs = sort(subs);

AllScene = table();
AllPairs = table();

perSub = table('Size',[0 8], 'VariableTypes', ...
    {'string','string','string','double','double','double','double','double'}, ...
    'VariableNames', {'subject_id','SportFreq','Experience','Order','n_pairs','mean_delta_O_alpha','mean_gap_sec','mean_gray_dur'});

for k = 1:numel(subs)
    sid = string(subs{k});

    % scene_level
    f_scene = fullfile(fp_out, char(sid), 'csv', sprintf('%s_scene_level.csv', sid));
    if exist(f_scene, 'file')
        try
            T = readtable(f_scene);
            T.subject_id = repmat(sid, height(T), 1);
            AllScene = [AllScene; T]; %#ok<AGROW>
        catch ME
            fprintf(2, '[WARN] Failed reading scene_level for %s: %s\n', sid, ME.message);
        end
    end

    % pairs_check
    f_pairs = fullfile(fp_out, char(sid), 'qc', sprintf('%s_pairs_check.csv', sid));
    if exist(f_pairs, 'file')
        try
            P = readtable(f_pairs);
            P.subject_id = repmat(sid, height(P), 1);
            AllPairs = [AllPairs; P]; %#ok<AGROW>

            % subject-level recovery metrics (include subject group labels when present)
            sf = ""; ex = ""; od = NaN;
            if ismember('SportFreq', P.Properties.VariableNames)
                sf = string(P.SportFreq(1));
            end
            if ismember('Experience', P.Properties.VariableNames)
                ex = string(P.Experience(1));
            end
            if ismember('Order', P.Properties.VariableNames)
                od = double(P.Order(1));
            end

            perSub = [perSub; {sid, sf, ex, od, height(P), mean(P.delta_O_alpha,'omitnan'), ...
                mean(P.gap_sec,'omitnan'), mean(P.gray_dur,'omitnan')}]; %#ok<AGROW>
        catch ME
            fprintf(2, '[WARN] Failed reading pairs_check for %s: %s\n', sid, ME.message);
        end
    end
end

out = struct();

% Write merged tables
if ~isempty(AllScene)
    out.all_subjects_scene_level = fullfile(fp_sum, 'all_subjects_scene_level.csv');
    writetable(AllScene, out.all_subjects_scene_level);
end

if ~isempty(AllPairs)
    out.all_subjects_pairs_check = fullfile(fp_sum, 'all_subjects_pairs_check.csv');
    writetable(AllPairs, out.all_subjects_pairs_check);
end

if height(perSub) > 0
    out.per_subject_recovery = fullfile(fp_sum, 'per_subject_recovery_metrics.csv');
    writetable(perSub, out.per_subject_recovery);
end

% QC filtering + reports (optional)
qcOut = struct();
AllScene_qc = AllScene;
AllPairs_qc = AllPairs;
try
    if isfield(cfg,'qc_apply') && logical(cfg.qc_apply)
        qcOut = pipeline.qc_filter_and_report(fp_out, fp_sum, cfg, AllScene, AllPairs);
        if isfield(qcOut,'AllScene_qc'); AllScene_qc = qcOut.AllScene_qc; end
        if isfield(qcOut,'AllPairs_qc'); AllPairs_qc = qcOut.AllPairs_qc; end

        % Print per-scene valid N and excluded subject names into log
        if isfield(qcOut,'qc_scene_valid_counts') && exist(qcOut.qc_scene_valid_counts,'file')
            try
                Tcnt = readtable(qcOut.qc_scene_valid_counts, 'TextType','string');
                fprintf('\n=== QC Scene Valid Counts (after exclusions) ===\n');
                for i=1:height(Tcnt)
                    sid = Tcnt.scene_id(i);
                    nt = Tcnt.n_total(i);
                    nv = Tcnt.n_valid(i);
                    excl = "";
                    if ismember('excluded_subjects', Tcnt.Properties.VariableNames)
                        excl = string(Tcnt.excluded_subjects(i));
                    end
                    label = "";
                    if ismember('scene_name', Tcnt.Properties.VariableNames)
                        label = string(Tcnt.scene_name(i));
                    end
                    if strlength(label)>0
                        fprintf('Scene %02d (%s): valid=%d/%d | excluded: %s\n', sid, label, nv, nt, excl);
                    else
                        fprintf('Scene %02d: valid=%d/%d | excluded: %s\n', sid, nv, nt, excl);
                    end
                end
            catch ME2
                fprintf(2, '[WARN] Failed printing qc_scene_valid_counts: %s\n', ME2.message);
            end
        end
    end
catch ME
    fprintf(2, '[WARN] qc_filter_and_report failed: %s\n', ME.message);
end

% Optional: group-level figures (PNG) under summary/fig/ (raw)
try
    pipeline.plot_group_summaries(AllScene, fp_sum, cfg);
catch ME
    fprintf(2, '[WARN] plot_group_summaries failed: %s\n', ME.message);
end

% Optional: group-level figures (PNG) under summary/fig/ (QC-filtered)
try
    if ~isequal(AllScene_qc, AllScene)
        cfg2 = cfg; cfg2.group_plot_metrics = getfield_def(cfg,'group_plot_metrics',{"O_alpha","O_theta","O_beta","F_theta"});
        pipeline.plot_group_summaries(AllScene_qc, fp_sum, cfg2);
    end
catch ME
    fprintf(2, '[WARN] plot_group_summaries(QC) failed: %s\n', ME.message);
end

% Optional: group-level recovery figures from pairs_check (PNG) under summary/fig/ (raw)
try
    pipeline.plot_group_recovery_summaries(AllPairs, fp_sum, cfg);
catch ME
    fprintf(2, '[WARN] plot_group_recovery_summaries failed: %s\n', ME.message);
end

% Optional: group-level recovery figures from pairs_check (PNG) under summary/fig/ (QC-filtered)
try
    if ~isequal(AllPairs_qc, AllPairs)
        pipeline.plot_group_recovery_summaries(AllPairs_qc, fp_sum, cfg);
    end
catch ME
    fprintf(2, '[WARN] plot_group_recovery_summaries(QC) failed: %s\n', ME.message);
end

% Optional: group-level topoplots under summary/fig/ (uses per-subject topo_long; QC by subject exclusion)
try
    if isfield(qcOut,'Qsub')
        cfg.qc_include_subjects = string(qcOut.Qsub.subject_id(~qcOut.Qsub.exclude_scenelevel));
    end
    pipeline.plot_group_topoplots(fp_out, fp_sum, cfg);
catch ME
    fprintf(2, '[WARN] plot_group_topoplots failed: %s\n', ME.message);
end

% Paper-ready multi-panel figures under summary/paper_fig/ (raw)
try
    pipeline.plot_paper_figures(fp_sum, cfg);
catch ME
    fprintf(2, '[WARN] plot_paper_figures failed: %s\n', ME.message);
end

% Paper-ready multi-panel figures under summary/paper_fig_qc/ (QC-filtered)
try
    if ~isequal(AllScene_qc, AllScene)
        fp_scene = fullfile(fp_sum,'all_subjects_scene_level_qc.csv');
        if exist(fp_scene,'file')
            % temporarily swap file name by copying
            fp_qc_dir = fullfile(fp_sum,'paper_fig_qc');
            if ~exist(fp_qc_dir,'dir'); mkdir(fp_qc_dir); end
            % use plot_paper_figures but point it to fp_sum; it reads standard filenames.
            % so we temporarily write QC tables as the standard names within fp_qc_dir.
            copyfile(fullfile(fp_sum,'all_subjects_scene_level_qc.csv'), fullfile(fp_qc_dir,'all_subjects_scene_level.csv'));
            if exist(fullfile(fp_sum,'all_subjects_pairs_check_qc.csv'),'file')
                copyfile(fullfile(fp_sum,'all_subjects_pairs_check_qc.csv'), fullfile(fp_qc_dir,'all_subjects_pairs_check.csv'));
            end
            pipeline.plot_paper_figures(fp_qc_dir, cfg);
        end
    end
catch ME
    fprintf(2, '[WARN] plot_paper_figures(QC) failed: %s\n', ME.message);
end

% Methods snapshot (markdown)
try
    pipeline.write_methods_snapshot(fp_sum, cfg);
catch
end

% Denominator sensitivity (1-30 vs 1-40)
try
    pipeline.compare_denominator_sensitivity(fp_sum, cfg);
catch ME
    fprintf(2, '[WARN] compare_denominator_sensitivity failed: %s\n', ME.message);
end

fprintf('Batch summaries written to: %s\n', fp_sum);

% Close summary-stage diary if we started it
try
    if didStartDiary
        fprintf('===== summarize_bandpower_outputs finished: %s =====\n', datestr(now,31));
        diary off;
        % Copy to pipeline_latest.log for convenience
        try
            copyfile(logPathRun, logPathLatest, 'f');
        catch
        end
    end
catch
end

end

function v = getfield_def(s, f, d)
try
    if isfield(s,f)
        v = s.(f);
    else
        v = d;
    end
catch
    v = d;
end
end
