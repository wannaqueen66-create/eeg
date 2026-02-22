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

% Optional: group-level figures (PNG) under summary/fig/
try
    pipeline.plot_group_summaries(AllScene, fp_sum, cfg);
catch ME
    fprintf(2, '[WARN] plot_group_summaries failed: %s\n', ME.message);
end

% Optional: group-level recovery figures from pairs_check (PNG) under summary/fig/
try
    pipeline.plot_group_recovery_summaries(AllPairs, fp_sum, cfg);
catch ME
    fprintf(2, '[WARN] plot_group_recovery_summaries failed: %s\n', ME.message);
end

% Optional: group-level topoplots under summary/fig/
try
    pipeline.plot_group_topoplots(fp_out, fp_sum, cfg);
catch ME
    fprintf(2, '[WARN] plot_group_topoplots failed: %s\n', ME.message);
end

fprintf('Batch summaries written to: %s\n', fp_sum);

end
