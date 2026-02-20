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

perSub = table('Size',[0 6], 'VariableTypes', ...
    {'string','double','double','double','double','double'}, ...
    'VariableNames', {'subject_id','n_pairs','mean_delta_O_alpha','mean_gap_sec','mean_view_dur','mean_gray_dur'});

for k = 1:numel(subs)
    sid = string(subs{k});

    % scene_level
    f_scene = fullfile(fp_out, char(sid), 'csv', sprintf('%s_scene_level.csv', sid));
    if exist(f_scene, 'file')
        try
            T = readtable(f_scene);
            T.subject_id = repmat(sid, height(T), 1);
            AllScene = [AllScene; T]; %#ok<AGROW>
        catch
        end
    end

    % pairs_check
    f_pairs = fullfile(fp_out, char(sid), 'qc', sprintf('%s_pairs_check.csv', sid));
    if exist(f_pairs, 'file')
        try
            P = readtable(f_pairs);
            P.subject_id = repmat(sid, height(P), 1);
            AllPairs = [AllPairs; P]; %#ok<AGROW>

            perSub = [perSub; {sid, height(P), mean(P.delta_O_alpha,'omitnan'), ...
                mean(P.gap_sec,'omitnan'), mean(P.view_dur,'omitnan'), mean(P.gray_dur,'omitnan')}]; %#ok<AGROW>
        catch
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

fprintf('Batch summaries written to: %s\n', fp_sum);

end
