function out = analyze_eye_eeg_links(AllSceneEye, fp_sum, cfg, tag)
%ANALYZE_EYE_EEG_LINKS Minimal EEG x eye exploratory analysis.
%
% Outputs under:
%   <summary>/analysis-eye/<tag>/
%     eye_eeg_correlation_screen.csv
%     eye_eeg_model_ready.csv
%     eye_eeg_lmm_fixed_effects_<metric>.csv
%     eye_eeg_lmm_anova_<metric>.csv
%     eye_eeg_report.md

if nargin < 4 || isempty(tag)
    tag = 'raw';
end
if nargin < 3
    cfg = struct();
end
out = struct();
if isempty(AllSceneEye) || ~istable(AllSceneEye)
    return;
end

req = {'subject_id','WWR','Complexity'};
for i = 1:numel(req)
    if ~ismember(req{i}, AllSceneEye.Properties.VariableNames)
        warning('analyze_eye_eeg_links: missing required column %s. Skipping.', req{i});
        return;
    end
end

fp_root = fullfile(fp_sum, 'analysis-eye', char(tag));
fp_tbl = fullfile(fp_root, 'tables');
fp_rep = fullfile(fp_root, 'reports');
if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end

% Default metric sets
EEG = {'F_theta','O_alpha','O_beta'};
EYE = {'eye_tracking_ratio','eye_view_blink_rate_per_min','eye_view_sacc_rate_per_min','eye_mean_pupil_mm'};
try
    if isfield(cfg,'eye_analysis_eeg_metrics') && ~isempty(cfg.eye_analysis_eeg_metrics)
        EEG = cellstr(string(cfg.eye_analysis_eeg_metrics));
    end
catch
end
try
    if isfield(cfg,'eye_analysis_eye_metrics') && ~isempty(cfg.eye_analysis_eye_metrics)
        EYE = cellstr(string(cfg.eye_analysis_eye_metrics));
    end
catch
end
EEG = EEG(ismember(EEG, AllSceneEye.Properties.VariableNames));
EYE = EYE(ismember(EYE, AllSceneEye.Properties.VariableNames));
if isempty(EEG) || isempty(EYE)
    warning('analyze_eye_eeg_links: no usable EEG or eye metrics found. Skipping.');
    return;
end

% Build ready table
T = AllSceneEye;
T.subject_id = string(T.subject_id);
T.Subject = categorical(strtrim(T.subject_id));
if ~iscategorical(T.Complexity); T.Complexity = categorical(double(T.Complexity)); end
if isnumeric(T.WWR)
    T.WWR = categorical(double(T.WWR));
else
    T.WWR = categorical(string(T.WWR));
end
out.model_ready_csv = fullfile(fp_tbl, 'eye_eeg_model_ready.csv');
writetable(T, out.model_ready_csv);

% Correlation screen
Rows = {};
for i = 1:numel(EEG)
    for j = 1:numel(EYE)
        xe = T.(EEG{i});
        ye = T.(EYE{j});
        try
            [r,p] = corr(double(xe), double(ye), 'Rows','complete', 'Type','Spearman');
        catch
            r = NaN; p = NaN;
        end
        n = sum(~isnan(double(xe)) & ~isnan(double(ye)));
        Rows(end+1,:) = {EEG{i}, EYE{j}, n, r, p}; %#ok<AGROW>
    end
end
C = cell2table(Rows, 'VariableNames', {'eeg_metric','eye_metric','n','spearman_r','p'});
out.correlation_csv = fullfile(fp_tbl, 'eye_eeg_correlation_screen.csv');
writetable(C, out.correlation_csv);

% Minimal LMMs: EEG ~ WWR*Complexity + eye_tracking_ratio + blink_rate + sacc_rate + (1|Subject)
if exist('fitlme','file') == 2
    lmmRows = {};
    for i = 1:numel(EEG)
        dv = EEG{i};
        varsNeed = [{'Subject','WWR','Complexity',dv}, intersect({'eye_tracking_ratio','eye_view_blink_rate_per_min','eye_view_sacc_rate_per_min','eye_mean_pupil_mm'}, T.Properties.VariableNames, 'stable')];
        K = T(:, varsNeed);
        K = rmmissing(K);
        if height(K) < 12
            continue;
        end
        f = sprintf('%s ~ WWR*Complexity', dv);
        if ismember('eye_tracking_ratio', K.Properties.VariableNames)
            f = [f ' + eye_tracking_ratio'];
        end
        if ismember('eye_view_blink_rate_per_min', K.Properties.VariableNames)
            f = [f ' + eye_view_blink_rate_per_min'];
        end
        if ismember('eye_view_sacc_rate_per_min', K.Properties.VariableNames)
            f = [f ' + eye_view_sacc_rate_per_min'];
        end
        if ismember('eye_mean_pupil_mm', K.Properties.VariableNames)
            f = [f ' + eye_mean_pupil_mm'];
        end
        f = [f ' + (1|Subject)'];
        try
            lme = fitlme(K, f);
            FE = lme.Coefficients;
            AN = anova(lme);
            fp_fe = fullfile(fp_tbl, sprintf('eye_eeg_lmm_fixed_effects_%s.csv', dv));
            fp_an = fullfile(fp_tbl, sprintf('eye_eeg_lmm_anova_%s.csv', dv));
            writetable(FE, fp_fe);
            writetable(AN, fp_an);
            lmmRows(end+1,:) = {dv, height(K), f, fp_fe, fp_an}; %#ok<AGROW>
        catch ME
            warning('analyze_eye_eeg_links: fitlme failed for %s: %s', dv, ME.message);
        end
    end
    if ~isempty(lmmRows)
        L = cell2table(lmmRows, 'VariableNames', {'eeg_metric','n_rows','formula','fixed_effects_csv','anova_csv'});
        out.lmm_index_csv = fullfile(fp_tbl, 'eye_eeg_lmm_index.csv');
        writetable(L, out.lmm_index_csv);
    end
end

% Markdown report
fp_md = fullfile(fp_rep, 'eye_eeg_report.md');
fid = fopen(fp_md,'w');
if fid ~= -1
    fprintf(fid, '# EEG × Eye Exploratory Report\n\n');
    fprintf(fid, '- Generated: %s\n', datestr(now,31));
    fprintf(fid, '- Data tag: `%s`\n\n', tag);

    fprintf(fid, '## Metrics used\n\n');
    fprintf(fid, '- EEG metrics: %s\n', strjoin(string(EEG), ', '));
    fprintf(fid, '- Eye metrics: %s\n\n', strjoin(string(EYE), ', '));

    fprintf(fid, '## Correlation screen\n\n');
    fprintf(fid, '- CSV: `%s`\n\n', out.correlation_csv);
    C2 = sortrows(C, 'p', 'ascend');
    topN = min(10, height(C2));
    fprintf(fid, '| eeg_metric | eye_metric | n | spearman_r | p |\n');
    fprintf(fid, '|---|---|---:|---:|---:|\n');
    for i = 1:topN
        fprintf(fid, '| %s | %s | %d | %.3f | %.4g |\n', string(C2.eeg_metric(i)), string(C2.eye_metric(i)), C2.n(i), C2.spearman_r(i), C2.p(i));
    end
    fprintf(fid, '\n');

    fprintf(fid, '## Minimal modeling template used\n\n');
    fprintf(fid, '- baseline idea: EEG ~ WWR*Complexity + eye_tracking_ratio + eye_view_blink_rate_per_min + eye_view_sacc_rate_per_min + (1|Subject)\n');
    if isfield(out,'lmm_index_csv') && exist(out.lmm_index_csv,'file')
        fprintf(fid, '- LMM index: `%s`\n', out.lmm_index_csv);
    else
        fprintf(fid, '- LMM results: not generated (likely Stats toolbox unavailable or too few complete rows).\n');
    end
    fprintf(fid, '\n');

    fprintf(fid, '## Interpretation notes\n\n');
    fprintf(fid, '- Use blink/saccade/tracking metrics first as QC/support variables, not as direct neural interpretations.\n');
    fprintf(fid, '- Pupil is better treated as state/attention support than artifact evidence.\n');
    fprintf(fid, '- Treat these results as exploratory screening before formal multimodal claims.\n');
    fclose(fid);
end
out.report_md = fp_md;
end
