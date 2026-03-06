function out = write_eye_qc_report_md(fp_sum, AllSceneEye, cfg)
%WRITE_EYE_QC_REPORT_MD Write markdown + CSV summaries for eye-derived QC flags.

out = struct();
if nargin < 2 || isempty(AllSceneEye) || ~istable(AllSceneEye)
    return;
end

req = {'subject_id','scene_id'};
for i = 1:numel(req)
    if ~ismember(req{i}, AllSceneEye.Properties.VariableNames)
        return;
    end
end

fp_rep = fp_sum;
try
    if exist('pipeline.get_report_dir','file')==2
        fp_rep = pipeline.get_report_dir(fp_sum, cfg);
    end
catch
end
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end

fp_tbl = pipeline.get_table_dir(fp_sum, cfg, 'merged_raw');
if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end

T = AllSceneEye;
T.subject_id = string(T.subject_id);

flagVars = intersect({ ...
    'eye_qc_low_tracking_flag', 'eye_qc_low_validity_flag', 'eye_qc_low_openness_flag', ...
    'eye_qc_high_blink_rate_flag', 'eye_qc_high_blink_burden_flag', 'eye_qc_high_sacc_rate_flag', ...
    'eye_qc_needs_review', 'eye_qc_severe_flag', 'eye_qc_flag_count'}, T.Properties.VariableNames, 'stable');
if isempty(flagVars)
    return;
end

% Subject-level summary
[Gsub, sub] = findgroups(T.subject_id);
Tsub = table(sub, splitapply(@numel, T.subject_id, Gsub), 'VariableNames', {'subject_id','n_scenes'});
if ismember('eye_qc_needs_review', T.Properties.VariableNames)
    Tsub.n_needs_review = splitapply(@(x) sum(double(x)>0), T.eye_qc_needs_review, Gsub);
else
    Tsub.n_needs_review = zeros(height(Tsub),1);
end
if ismember('eye_qc_severe_flag', T.Properties.VariableNames)
    Tsub.n_severe = splitapply(@(x) sum(double(x)>0), T.eye_qc_severe_flag, Gsub);
else
    Tsub.n_severe = zeros(height(Tsub),1);
end
if ismember('eye_qc_flag_count', T.Properties.VariableNames)
    Tsub.mean_flag_count = splitapply(@(x) mean(double(x),'omitnan'), T.eye_qc_flag_count, Gsub);
else
    Tsub.mean_flag_count = nan(height(Tsub),1);
end
Tsub.frac_needs_review = Tsub.n_needs_review ./ max(Tsub.n_scenes,1);
Tsub.frac_severe = Tsub.n_severe ./ max(Tsub.n_scenes,1);
out.eye_qc_subject_summary_csv = fullfile(fp_tbl, 'eye_qc_subject_summary.csv');
writetable(Tsub, out.eye_qc_subject_summary_csv);

% Scene-level summary
[Gsc, sc] = findgroups(T.scene_id);
Tsc = table(sc, splitapply(@numel, T.scene_id, Gsc), 'VariableNames', {'scene_id','n_rows'});
if ismember('scene_name', T.Properties.VariableNames)
    Tsc.scene_name = splitapply(@(x) string(x(find(strlength(string(x))>0,1,'first'))), string(T.scene_name), Gsc);
end
if ismember('eye_qc_needs_review', T.Properties.VariableNames)
    Tsc.n_needs_review = splitapply(@(x) sum(double(x)>0), T.eye_qc_needs_review, Gsc);
else
    Tsc.n_needs_review = zeros(height(Tsc),1);
end
if ismember('eye_qc_severe_flag', T.Properties.VariableNames)
    Tsc.n_severe = splitapply(@(x) sum(double(x)>0), T.eye_qc_severe_flag, Gsc);
else
    Tsc.n_severe = zeros(height(Tsc),1);
end
if ismember('eye_qc_flag_count', T.Properties.VariableNames)
    Tsc.mean_flag_count = splitapply(@(x) mean(double(x),'omitnan'), T.eye_qc_flag_count, Gsc);
else
    Tsc.mean_flag_count = nan(height(Tsc),1);
end
Tsc.frac_needs_review = Tsc.n_needs_review ./ max(Tsc.n_rows,1);
Tsc.frac_severe = Tsc.n_severe ./ max(Tsc.n_rows,1);
out.eye_qc_scene_summary_csv = fullfile(fp_tbl, 'eye_qc_scene_summary.csv');
writetable(Tsc, out.eye_qc_scene_summary_csv);

% Human-readable report
fp_md = fullfile(fp_rep, 'eye_qc_report.md');
fid = fopen(fp_md,'w');
if fid ~= -1
    fprintf(fid, '# Eye QC Report (Auto-generated)\n\n');
    fprintf(fid, '- Generated: %s\n', datestr(now,31));
    fprintf(fid, '- Rows in merged scene table: %d\n', height(T));
    if ismember('eye_qc_needs_review', T.Properties.VariableNames)
        fprintf(fid, '- Scenes needing review: %d (%.1f%%)\n', sum(double(T.eye_qc_needs_review)>0), 100*mean(double(T.eye_qc_needs_review)>0,'omitnan'));
    end
    if ismember('eye_qc_severe_flag', T.Properties.VariableNames)
        fprintf(fid, '- Scenes with severe eye QC flags: %d (%.1f%%)\n\n', sum(double(T.eye_qc_severe_flag)>0), 100*mean(double(T.eye_qc_severe_flag)>0,'omitnan'));
    else
        fprintf(fid, '\n');
    end

    fprintf(fid, '## Flag columns included\n\n');
    for i = 1:numel(flagVars)
        fprintf(fid, '- `%s`\n', flagVars{i});
    end
    fprintf(fid, '\n');

    fprintf(fid, '## Per-scene summary\n\n');
    fprintf(fid, '| scene_id | scene_name | needs_review / total | severe / total | mean_flag_count |\n');
    fprintf(fid, '|---:|---|---:|---:|---:|\n');
    for i = 1:height(Tsc)
        sname = "";
        if ismember('scene_name', Tsc.Properties.VariableNames)
            sname = string(Tsc.scene_name(i));
        end
        fprintf(fid, '| %d | %s | %d / %d | %d / %d | %.2f |\n', ...
            Tsc.scene_id(i), sname, Tsc.n_needs_review(i), Tsc.n_rows(i), Tsc.n_severe(i), Tsc.n_rows(i), Tsc.mean_flag_count(i));
    end
    fprintf(fid, '\n');

    fprintf(fid, '## Most flagged subjects\n\n');
    Tsub2 = sortrows(Tsub, {'n_severe','n_needs_review','mean_flag_count'}, {'descend','descend','descend'});
    topN = min(15, height(Tsub2));
    fprintf(fid, '| subject_id | n_scenes | needs_review | severe | mean_flag_count |\n');
    fprintf(fid, '|---|---:|---:|---:|---:|\n');
    for i = 1:topN
        fprintf(fid, '| %s | %d | %d | %d | %.2f |\n', string(Tsub2.subject_id(i)), Tsub2.n_scenes(i), Tsub2.n_needs_review(i), Tsub2.n_severe(i), Tsub2.mean_flag_count(i));
    end
    fprintf(fid, '\n');

    fprintf(fid, '## Suggested first covariates\n\n');
    fprintf(fid, '- `eye_tracking_ratio`\n');
    fprintf(fid, '- `eye_view_blink_rate_per_min`\n');
    fprintf(fid, '- `eye_view_sacc_rate_per_min`\n');
    fprintf(fid, '- optional: `eye_view_blink_burden_pct`, `eye_mean_pupil_mm`\n\n');

    fprintf(fid, '## Notes\n\n');
    fprintf(fid, '- These flags are exploratory QC helpers, not hard exclusion rules.\n');
    fprintf(fid, '- Use them to prioritize manual review and interpret suspicious EEG patterns more cautiously.\n');
    fclose(fid);
end
out.eye_qc_report_md = fp_md;
end
