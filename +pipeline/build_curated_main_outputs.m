function out = build_curated_main_outputs(fp_sum, cfg, AllScene, AllPairs, AllScene_qc, AllPairs_qc)
%BUILD_CURATED_MAIN_OUTPUTS Build a cleaner top-level output surface for the redesigned main branch.
%
% This function does NOT replace existing task outputs.
% It creates a curated, easier-to-read layer organized as:
%   <batch>/descriptive/{overall,experience}/...
%   <batch>/inferential/{overall,experience}/...
%
% The goal is to give the repository a clearer public-facing structure while
% preserving the current detailed task outputs underneath.
%
% Current redesign priority:
% - overall
% - experience
%
% SportFreq remains available in the detailed task outputs, but is not exposed
% as a first-class curated surface in the redesigned main branch.

if nargin < 6; AllPairs_qc = table(); end
if nargin < 5; AllScene_qc = table(); end
if nargin < 4; AllPairs = table(); end
if nargin < 3; AllScene = table(); end
if nargin < 2; cfg = struct(); end

out = struct();

if isempty(fp_sum) || ~exist(fp_sum,'dir')
    return;
end

% Expect fp_sum to be the batch root under the staged layout.
fp_desc = fullfile(fp_sum, 'descriptive');
fp_inf  = fullfile(fp_sum, 'inferential');
mkdir_if_needed(fp_desc);
mkdir_if_needed(fp_inf);

% Batch-level guide for the redesigned main output surface
try
    fp_batch_rep = fullfile(fp_sum, 'reports');
    mkdir_if_needed(fp_batch_rep);
    fid = fopen(fullfile(fp_batch_rep, 'README_CURATED_MAIN.md'),'w');
    if fid>0
        fprintf(fid, '# Curated Main Output Surface\n\n');
        fprintf(fid, 'Priority reading order for the redesigned main branch:\n');
        fprintf(fid, '1. `descriptive/overall/`\n');
        fprintf(fid, '2. `descriptive/experience/`\n');
        fprintf(fid, '3. `inferential/overall/`\n');
        fprintf(fid, '4. `inferential/experience/`\n\n');
        fprintf(fid, 'The detailed task-based outputs remain available under `analysis/`, but the redesigned main branch intentionally prioritizes only two visible branches by default: `overall` and `experience`.\n\n');
        fprintf(fid, 'Note: SportFreq analyses are still present in detailed outputs, but are currently treated as secondary rather than first-class curated outputs.\n');
        fclose(fid);
    end
catch
end
try
    pipeline.write_curated_readme_index(fp_sum);
catch
end

% ---------- descriptive / overall ----------
fp_do = fullfile(fp_desc, 'overall');
fp_do_tbl = fullfile(fp_do, 'tables');
fp_do_fig = fullfile(fp_do, 'figures');
fp_do_rep = fullfile(fp_do, 'report');
mkdir_if_needed(fp_do_tbl); mkdir_if_needed(fp_do_fig); mkdir_if_needed(fp_do_rep);
out.descriptive_overall = fp_do;

try
    if ~isempty(AllScene)
        writetable(AllScene, fullfile(fp_do_tbl, 'scene_level_overall_raw.csv'));
        Tsum = summarize_overall_scene(AllScene, cfg);
        if ~isempty(Tsum)
            writetable(Tsum, fullfile(fp_do_tbl, 'overall_scene_metric_means_raw.csv'));
            plot_overall_metric_bar(Tsum, 'raw', fp_do_fig, cfg);
        end
        Tcond = summarize_overall_by_factors(AllScene, cfg);
        if ~isempty(Tcond)
            writetable(Tcond, fullfile(fp_do_tbl, 'overall_scene_metric_means_by_WWR_Complexity_raw.csv'));
            plot_overall_factor_grid(Tcond, 'raw', fp_do_fig, cfg);
        end
        Ttrial = summarize_trialindex_overall(AllScene, cfg);
        if ~isempty(Ttrial)
            writetable(Ttrial, fullfile(fp_do_tbl, 'overall_trialindex_neural_response_raw.csv'));
            plot_trialindex_descriptive_overall(Ttrial, 'raw', fp_do_fig, cfg);
        end
        TtrialOrd = summarize_trialindex_by_order(AllScene, cfg);
        if ~isempty(TtrialOrd)
            writetable(TtrialOrd, fullfile(fp_do_tbl, 'overall_trialindex_neural_response_by_order_raw.csv'));
            plot_trialindex_descriptive_by_order(TtrialOrd, 'raw', fp_do_fig, cfg);
        end
    end
    if ~isempty(AllScene_qc)
        writetable(AllScene_qc, fullfile(fp_do_tbl, 'scene_level_overall_qc.csv'));
        Tsum = summarize_overall_scene(AllScene_qc, cfg);
        if ~isempty(Tsum)
            writetable(Tsum, fullfile(fp_do_tbl, 'overall_scene_metric_means_qc.csv'));
            plot_overall_metric_bar(Tsum, 'qc', fp_do_fig, cfg);
        end
        Tcond = summarize_overall_by_factors(AllScene_qc, cfg);
        if ~isempty(Tcond)
            writetable(Tcond, fullfile(fp_do_tbl, 'overall_scene_metric_means_by_WWR_Complexity_qc.csv'));
            plot_overall_factor_grid(Tcond, 'qc', fp_do_fig, cfg);
        end
        Ttrial = summarize_trialindex_overall(AllScene_qc, cfg);
        if ~isempty(Ttrial)
            writetable(Ttrial, fullfile(fp_do_tbl, 'overall_trialindex_neural_response_qc.csv'));
            plot_trialindex_descriptive_overall(Ttrial, 'qc', fp_do_fig, cfg);
        end
        TtrialOrd = summarize_trialindex_by_order(AllScene_qc, cfg);
        if ~isempty(TtrialOrd)
            writetable(TtrialOrd, fullfile(fp_do_tbl, 'overall_trialindex_neural_response_by_order_qc.csv'));
            plot_trialindex_descriptive_by_order(TtrialOrd, 'qc', fp_do_fig, cfg);
        end
    end
    if ~isempty(AllPairs)
        writetable(AllPairs, fullfile(fp_do_tbl, 'pairs_overall_raw.csv'));
        ToverRec = summarize_recovery_overall(AllPairs);
        if ~isempty(ToverRec)
            writetable(ToverRec, fullfile(fp_do_tbl, 'overall_recovery_means_raw.csv'));
            plot_recovery_overall(ToverRec, 'raw', fp_do_fig, cfg);
        end
    end
    if ~isempty(AllPairs_qc)
        writetable(AllPairs_qc, fullfile(fp_do_tbl, 'pairs_overall_qc.csv'));
        ToverRec = summarize_recovery_overall(AllPairs_qc);
        if ~isempty(ToverRec)
            writetable(ToverRec, fullfile(fp_do_tbl, 'overall_recovery_means_qc.csv'));
            plot_recovery_overall(ToverRec, 'qc', fp_do_fig, cfg);
        end
    end
catch
end

try
    fid = fopen(fullfile(fp_do_rep, 'README.md'),'w');
    if fid>0
        fprintf(fid, '# Descriptive / Overall\n\n');
        fprintf(fid, '- `scene_level_overall_raw.csv`: full-sample scene-level descriptive table\n');
        fprintf(fid, '- `scene_level_overall_qc.csv`: QC-filtered scene-level descriptive table\n');
        fprintf(fid, '- `overall_scene_metric_means_raw.csv`: core metric means across the full sample\n');
        fprintf(fid, '- `overall_scene_metric_means_qc.csv`: QC-filtered core metric means across the full sample\n');
        fprintf(fid, '- `overall_scene_metric_means_by_WWR_Complexity_raw.csv`: descriptive means by WWR × Complexity\n');
        fprintf(fid, '- `overall_scene_metric_means_by_WWR_Complexity_qc.csv`: QC-filtered descriptive means by WWR × Complexity\n');
        fprintf(fid, '- `overall_trialindex_neural_response_*.csv`: trial-by-trial neural response summaries across TrialIndex 1–12\n');
        fprintf(fid, '- `overall_trialindex_neural_response_by_order_*.csv`: trial-by-trial neural response summaries split by counterbalanced Order\n');
        fprintf(fid, '- `pairs_overall_raw.csv`: full-sample recovery/pair descriptive table\n');
        fprintf(fid, '- `pairs_overall_qc.csv`: QC-filtered recovery/pair descriptive table\n');
        fprintf(fid, '- `overall_recovery_means_*.csv`: overall recovery summary table\n');
        fprintf(fid, '- `overall_metric_bar_*.png`: overall core metric bar figure\n');
        fprintf(fid, '- `overall_recovery_bar_*.png`: overall recovery bar figure\n');
        fprintf(fid, '- `overall_factor_grid_*.png`: overall WWR × Complexity descriptive grid\n');
        fprintf(fid, '- `overall_trialindex_response_*.png`: overall trial-by-trial neural response figure\n');
        fprintf(fid, '- `overall_trialindex_response_by_order_*.png`: overall trial-by-trial neural response figure split by Order\n');
        fprintf(fid, '\nThis folder is intended as the simplest entry point for reading descriptive EEG results at the full-sample level.\n');
        fclose(fid);
    end
catch
end

% ---------- descriptive / experience ----------
fp_de = fullfile(fp_desc, 'experience');
fp_de_tbl = fullfile(fp_de, 'tables');
fp_de_fig = fullfile(fp_de, 'figures');
fp_de_rep = fullfile(fp_de, 'report');
mkdir_if_needed(fp_de_tbl); mkdir_if_needed(fp_de_fig); mkdir_if_needed(fp_de_rep);
out.descriptive_experience = fp_de;

expColScene = pick_group_col(AllScene, 'ExperienceGroup', 'Experience');
expColPairs = pick_group_col(AllPairs, 'ExperienceGroup', 'Experience');

try
    if ~isempty(AllScene) && strlength(expColScene)>0
        S = AllScene;
        G = normalize_high_low_local(S.(char(expColScene)));
        S.ExperienceGroupCurated = G;
        S = S(strlength(G)>0,:);
        if ~isempty(S)
            writetable(S, fullfile(fp_de_tbl, 'scene_level_experience_raw.csv'));
            Tsum = summarize_scene_by_experience(S, cfg);
            if ~isempty(Tsum)
                writetable(Tsum, fullfile(fp_de_tbl, 'experience_scene_metric_means_raw.csv'));
                plot_experience_metric_bar(Tsum, 'raw', fp_de_fig, cfg);
            end
            Tcond = summarize_scene_by_experience_factors(S, cfg);
            if ~isempty(Tcond)
                writetable(Tcond, fullfile(fp_de_tbl, 'experience_scene_metric_means_by_WWR_Complexity_raw.csv'));
                plot_experience_factor_grid(Tcond, 'raw', fp_de_fig, cfg);
            end
            Ttrial = summarize_trialindex_experience(S, cfg);
            if ~isempty(Ttrial)
                writetable(Ttrial, fullfile(fp_de_tbl, 'experience_trialindex_neural_response_raw.csv'));
                plot_trialindex_descriptive_experience(Ttrial, 'raw', fp_de_fig, cfg);
            end
        end
    end
    if ~isempty(AllScene_qc)
        expColSceneQc = pick_group_col(AllScene_qc, 'ExperienceGroup', 'Experience');
        if strlength(expColSceneQc)>0
            S = AllScene_qc;
            G = normalize_high_low_local(S.(char(expColSceneQc)));
            S.ExperienceGroupCurated = G;
            S = S(strlength(G)>0,:);
            if ~isempty(S)
                writetable(S, fullfile(fp_de_tbl, 'scene_level_experience_qc.csv'));
                Tsum = summarize_scene_by_experience(S, cfg);
                if ~isempty(Tsum)
                    writetable(Tsum, fullfile(fp_de_tbl, 'experience_scene_metric_means_qc.csv'));
                    plot_experience_metric_bar(Tsum, 'qc', fp_de_fig, cfg);
                end
                Tcond = summarize_scene_by_experience_factors(S, cfg);
                if ~isempty(Tcond)
                    writetable(Tcond, fullfile(fp_de_tbl, 'experience_scene_metric_means_by_WWR_Complexity_qc.csv'));
                    plot_experience_factor_grid(Tcond, 'qc', fp_de_fig, cfg);
                end
                Ttrial = summarize_trialindex_experience(S, cfg);
                if ~isempty(Ttrial)
                    writetable(Ttrial, fullfile(fp_de_tbl, 'experience_trialindex_neural_response_qc.csv'));
                    plot_trialindex_descriptive_experience(Ttrial, 'qc', fp_de_fig, cfg);
                end
            end
        end
    end
    if ~isempty(AllPairs) && strlength(expColPairs)>0
        P = AllPairs;
        G = normalize_high_low_local(P.(char(expColPairs)));
        P.ExperienceGroupCurated = G;
        P = P(strlength(G)>0,:);
        if ~isempty(P)
            writetable(P, fullfile(fp_de_tbl, 'pairs_experience_raw.csv'));
            Trec = summarize_recovery_by_experience(P);
            if ~isempty(Trec)
                writetable(Trec, fullfile(fp_de_tbl, 'experience_recovery_means_raw.csv'));
                plot_recovery_experience(Trec, 'raw', fp_de_fig, cfg);
            end
        end
    end
    if ~isempty(AllPairs_qc)
        expColPairsQc = pick_group_col(AllPairs_qc, 'ExperienceGroup', 'Experience');
        if strlength(expColPairsQc)>0
            P = AllPairs_qc;
            G = normalize_high_low_local(P.(char(expColPairsQc)));
            P.ExperienceGroupCurated = G;
            P = P(strlength(G)>0,:);
            if ~isempty(P)
                writetable(P, fullfile(fp_de_tbl, 'pairs_experience_qc.csv'));
                Trec = summarize_recovery_by_experience(P);
                if ~isempty(Trec)
                    writetable(Trec, fullfile(fp_de_tbl, 'experience_recovery_means_qc.csv'));
                    plot_recovery_experience(Trec, 'qc', fp_de_fig, cfg);
                end
            end
        end
    end
catch
end

% Copy a few already-generated experience-oriented figures if present
try
    copy_matching_figs(fullfile(fp_sum,'figures','group','raw'), fp_de_fig, 'experience');
    copy_matching_figs(fullfile(fp_sum,'figures','group','qc'),  fp_de_fig, 'experience');
    copy_matching_figs(fullfile(fp_sum,'figures','recovery','raw'), fp_de_fig, 'experience');
    copy_matching_figs(fullfile(fp_sum,'figures','recovery','qc'),  fp_de_fig, 'experience');
catch
end

try
    fid = fopen(fullfile(fp_de_rep, 'README.md'),'w');
    if fid>0
        fprintf(fid, '# Descriptive / Experience\n\n');
        fprintf(fid, '- `scene_level_experience_*.csv`: scene-level descriptive tables with valid experience labels\n');
        fprintf(fid, '- `experience_scene_metric_means_*.csv`: metric summaries by Experience group (and Complexity when available)\n');
        fprintf(fid, '- `experience_scene_metric_means_by_WWR_Complexity_*.csv`: metric summaries by Experience × WWR × Complexity\n');
        fprintf(fid, '- `experience_trialindex_neural_response_*.csv`: trial-by-trial neural response summaries by Experience group\n');
        fprintf(fid, '- `pairs_experience_*.csv`: recovery tables with valid experience labels\n');
        fprintf(fid, '- `experience_recovery_means_*.csv`: recovery summaries by Experience group\n');
        fprintf(fid, '- `experience_metric_bar_*.png`: experience-group core metric bar figure\n');
        fprintf(fid, '- `experience_recovery_bar_*.png`: experience-group recovery bar figure\n');
        fprintf(fid, '- `experience_factor_grid_*.png`: experience-group WWR × Complexity descriptive grid\n');
        fprintf(fid, '- `experience_trialindex_response_*.png`: experience-group trial-by-trial neural response figure\n');
        fprintf(fid, '\nThis folder is intended as the default grouped descriptive layer for the redesigned main branch.\n');
        fclose(fid);
    end
catch
end

% ---------- inferential / overall ----------
fp_io = fullfile(fp_inf, 'overall');
fp_io_tbl = fullfile(fp_io, 'tables');
fp_io_fig = fullfile(fp_io, 'figures');
fp_io_rep = fullfile(fp_io, 'report');
mkdir_if_needed(fp_io_tbl); mkdir_if_needed(fp_io_fig); mkdir_if_needed(fp_io_rep);
out.inferential_overall = fp_io;

try
    build_inferential_overall_summary(fp_sum, fp_io_tbl, fp_io_fig, fp_io_rep);
catch
end
try
    if ~isempty(AllScene)
        build_inferential_overall_models(AllScene, cfg, 'raw', fp_io_tbl, fp_io_fig, fp_io_rep);
        build_inferential_overall_trend(AllScene, cfg, 'raw', fp_io_tbl, fp_io_fig, fp_io_rep);
    end
    if ~isempty(AllScene_qc)
        build_inferential_overall_models(AllScene_qc, cfg, 'qc', fp_io_tbl, fp_io_fig, fp_io_rep);
        build_inferential_overall_trend(AllScene_qc, cfg, 'qc', fp_io_tbl, fp_io_fig, fp_io_rep);
    end
    if ~isempty(AllPairs)
        build_inferential_overall_recovery(AllPairs, 'raw', fp_io_tbl, fp_io_fig, fp_io_rep);
    end
    if ~isempty(AllPairs_qc)
        build_inferential_overall_recovery(AllPairs_qc, 'qc', fp_io_tbl, fp_io_fig, fp_io_rep);
    end
catch
end

% ---------- inferential / experience ----------
fp_ie = fullfile(fp_inf, 'experience');
fp_ie_tbl = fullfile(fp_ie, 'tables');
fp_ie_fig = fullfile(fp_ie, 'figures');
fp_ie_rep = fullfile(fp_ie, 'report');
mkdir_if_needed(fp_ie_tbl); mkdir_if_needed(fp_ie_fig); mkdir_if_needed(fp_ie_rep);
out.inferential_experience = fp_ie;

% Mirror selected task-based experience outputs into a cleaner experience-facing surface.
try
    mirror_task_branch(fullfile(fp_sum,'analysis','task3_trialindex_lmm'), 'experience', fp_ie);
    mirror_task_branch(fullfile(fp_sum,'analysis','task4_core_lmm_suite'), 'experience', fp_ie);
    mirror_task_branch(fullfile(fp_sum,'analysis','task5_peakindex_invertedu'), 'experience', fp_ie);
    mirror_task_branch(fullfile(fp_sum,'analysis','task6_coremetric_special'), 'experience', fp_ie);
    mirror_task_branch(fullfile(fp_sum,'analysis','task7_individual_checks'), 'experience', fp_ie);
catch
end

try
    build_inferential_experience_summary(fp_sum, fp_ie_tbl, fp_ie_fig, fp_ie_rep);
catch
end
try
    if ~isempty(AllScene)
        build_inferential_experience_models(AllScene, cfg, 'raw', fp_ie_tbl, fp_ie_fig, fp_ie_rep);
        build_inferential_experience_trend(AllScene, cfg, 'raw', fp_ie_tbl, fp_ie_fig, fp_ie_rep);
    end
    if ~isempty(AllScene_qc)
        build_inferential_experience_models(AllScene_qc, cfg, 'qc', fp_ie_tbl, fp_ie_fig, fp_ie_rep);
        build_inferential_experience_trend(AllScene_qc, cfg, 'qc', fp_ie_tbl, fp_ie_fig, fp_ie_rep);
    end
    if ~isempty(AllPairs)
        build_inferential_experience_recovery(AllPairs, 'raw', fp_ie_tbl, fp_ie_fig, fp_ie_rep);
    end
    if ~isempty(AllPairs_qc)
        build_inferential_experience_recovery(AllPairs_qc, 'qc', fp_ie_tbl, fp_ie_fig, fp_ie_rep);
    end
catch
end

try
    write_curated_node_report(fp_do, 'descriptive', 'overall');
    write_curated_node_report(fp_de, 'descriptive', 'experience');
    write_curated_node_report(fp_io, 'inferential', 'overall');
    write_curated_node_report(fp_ie, 'inferential', 'experience');
catch
end

end

function mkdir_if_needed(p)
if ~exist(p,'dir'); mkdir(p); end
end

function Tsum = summarize_trialindex_overall(S, cfg)
Tsum = table();
if isempty(S) || ~ismember('subject_id', S.Properties.VariableNames)
    return;
end
trial = compute_trialindex_local(S);
if all(~isfinite(trial))
    return;
end
metrics = get_curated_metrics_local(S, cfg);
rows = {};
for mi = 1:numel(metrics)
    m = string(metrics{mi});
    if ~ismember(m, S.Properties.VariableNames), continue; end
    y = double(S.(char(m)));
    for ti = 1:12
        use = isfinite(trial) & trial==ti & isfinite(y);
        if ~any(use), continue; end
        subj = string(S.subject_id(use));
        y_use = y(use);
        subj_u = unique(subj);
        subj_means = nan(numel(subj_u),1);
        for si = 1:numel(subj_u)
            subj_means(si) = mean(y_use(subj==subj_u(si)), 'omitnan');
        end
        rows(end+1,:) = {m, ti, sum(use), numel(subj_u), mean(subj_means,'omitnan'), std(subj_means,'omitnan'), std(subj_means,'omitnan')/max(1,sqrt(numel(subj_u)))}; %#ok<AGROW>
    end
end
if isempty(rows), return; end
Tsum = cell2table(rows, 'VariableNames', {'metric','TrialIndex','n_rows','n_subjects','mean_value','sd_subject_mean','sem_subject_mean'});
end

function Tsum = summarize_trialindex_by_order(S, cfg)
Tsum = table();
if isempty(S) || ~ismember('Order', S.Properties.VariableNames) || ~ismember('subject_id', S.Properties.VariableNames)
    return;
end
trial = compute_trialindex_local(S);
if all(~isfinite(trial))
    return;
end
Ord = normalize_order_local(S.Order);
metrics = get_curated_metrics_local(S, cfg);
rows = {};
for mi = 1:numel(metrics)
    m = string(metrics{mi});
    if ~ismember(m, S.Properties.VariableNames), continue; end
    y = double(S.(char(m)));
    for ord = unique(Ord(:))'
        if strlength(ord)==0, continue; end
        for ti = 1:12
            use = isfinite(trial) & trial==ti & isfinite(y) & Ord==ord;
            if ~any(use), continue; end
            subj = string(S.subject_id(use));
            y_use = y(use);
            subj_u = unique(subj);
            subj_means = nan(numel(subj_u),1);
            for si = 1:numel(subj_u)
                subj_means(si) = mean(y_use(subj==subj_u(si)), 'omitnan');
            end
            rows(end+1,:) = {m, ord, ti, sum(use), numel(subj_u), mean(subj_means,'omitnan'), std(subj_means,'omitnan'), std(subj_means,'omitnan')/max(1,sqrt(numel(subj_u)))}; %#ok<AGROW>
        end
    end
end
if isempty(rows), return; end
Tsum = cell2table(rows, 'VariableNames', {'metric','Order','TrialIndex','n_rows','n_subjects','mean_value','sd_subject_mean','sem_subject_mean'});
end

function Tsum = summarize_trialindex_experience(S, cfg)
Tsum = table();
expCol = pick_group_col(S, 'ExperienceGroup', 'Experience');
if isempty(S) || strlength(expCol)==0 || ~ismember('subject_id', S.Properties.VariableNames)
    return;
end
trial = compute_trialindex_local(S);
if all(~isfinite(trial))
    return;
end
G = normalize_high_low_local(S.(char(expCol)));
metrics = get_curated_metrics_local(S, cfg);
rows = {};
for mi = 1:numel(metrics)
    m = string(metrics{mi});
    if ~ismember(m, S.Properties.VariableNames), continue; end
    y = double(S.(char(m)));
    for grp = ["Low","High"]
        for ti = 1:12
            use = isfinite(trial) & trial==ti & isfinite(y) & G==grp;
            if ~any(use), continue; end
            subj = string(S.subject_id(use));
            y_use = y(use);
            subj_u = unique(subj);
            subj_means = nan(numel(subj_u),1);
            for si = 1:numel(subj_u)
                subj_means(si) = mean(y_use(subj==subj_u(si)), 'omitnan');
            end
            rows(end+1,:) = {m, grp, ti, sum(use), numel(subj_u), mean(subj_means,'omitnan'), std(subj_means,'omitnan'), std(subj_means,'omitnan')/max(1,sqrt(numel(subj_u)))}; %#ok<AGROW>
        end
    end
end
if isempty(rows), return; end
Tsum = cell2table(rows, 'VariableNames', {'metric','ExperienceGroup','TrialIndex','n_rows','n_subjects','mean_value','sd_subject_mean','sem_subject_mean'});
end

function plot_trialindex_descriptive_overall(Tsum, tag, fp_fig, cfg)
if isempty(Tsum) || height(Tsum)==0, return; end
metrics = unique(string(Tsum.metric), 'stable');
if isempty(metrics), return; end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 1100 700]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
for i=1:min(4,numel(metrics))
    ax = nexttile(tl); hold(ax,'on'); style_axes(ax);
    M = Tsum(string(Tsum.metric)==metrics(i),:);
    errorbar(ax, double(M.TrialIndex), double(M.mean_value), double(M.sem_subject_mean), '-o', 'LineWidth',1.5, 'MarkerSize',4, 'Color',[0.15 0.35 0.7]);
    xlabel(ax,'TrialIndex'); ylabel(ax,'Mean neural response');
    title(ax, char(metrics(i)), 'Interpreter','none', 'FontWeight','normal');
    xlim(ax,[1 12]);
end
title(tl, sprintf('Overall trial-by-trial neural response [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('overall_trialindex_response_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function plot_trialindex_descriptive_experience(Tsum, tag, fp_fig, cfg)
if isempty(Tsum) || height(Tsum)==0, return; end
metrics = unique(string(Tsum.metric), 'stable');
if isempty(metrics), return; end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 1100 700]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
for i=1:min(4,numel(metrics))
    ax = nexttile(tl); hold(ax,'on'); style_axes(ax);
    M = Tsum(string(Tsum.metric)==metrics(i),:);
    for grp = ["Low","High"]
        Mg = M(string(M.ExperienceGroup)==grp,:);
        if isempty(Mg), continue; end
        if grp=="Low"
            cc = [0.18 0.49 0.72];
        else
            cc = [0.88 0.47 0.18];
        end
        errorbar(ax, double(Mg.TrialIndex), double(Mg.mean_value), double(Mg.sem_subject_mean), '-o', 'LineWidth',1.4, 'MarkerSize',4, 'Color',cc, 'DisplayName',char(grp));
    end
    xlabel(ax,'TrialIndex'); ylabel(ax,'Mean neural response');
    title(ax, char(metrics(i)), 'Interpreter','none', 'FontWeight','normal');
    xlim(ax,[1 12]); legend(ax,'Location','best');
end
title(tl, sprintf('Experience-group trial-by-trial neural response [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('experience_trialindex_response_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function plot_trialindex_descriptive_by_order(Tsum, tag, fp_fig, cfg)
if isempty(Tsum) || height(Tsum)==0, return; end
metrics = unique(string(Tsum.metric), 'stable');
if isempty(metrics), return; end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 1100 700]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
for i=1:min(4,numel(metrics))
    ax = nexttile(tl); hold(ax,'on'); style_axes(ax);
    M = Tsum(string(Tsum.metric)==metrics(i),:);
    ords = unique(string(M.Order), 'stable');
    cmap = lines(max(2,numel(ords)));
    for oi = 1:numel(ords)
        Mo = M(string(M.Order)==ords(oi),:);
        if isempty(Mo), continue; end
        errorbar(ax, double(Mo.TrialIndex), double(Mo.mean_value), double(Mo.sem_subject_mean), '-o', 'LineWidth',1.4, 'MarkerSize',4, 'Color',cmap(oi,:), 'DisplayName',char(ords(oi)));
    end
    xlabel(ax,'TrialIndex'); ylabel(ax,'Mean neural response');
    title(ax, char(metrics(i)), 'Interpreter','none', 'FontWeight','normal');
    xlim(ax,[1 12]); legend(ax,'Location','best');
end
title(tl, sprintf('Overall trial-by-trial neural response by Order [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('overall_trialindex_response_by_order_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function ord = normalize_order_local(x)
ord = repmat("", numel(x), 1);
s = strtrim(string(x));
for i=1:numel(s)
    if ~isnan(str2double(s(i)))
        ord(i) = "Order" + string(str2double(s(i)));
    elseif strlength(s(i))>0
        tok = regexp(char(s(i)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            ord(i) = "Order" + string(str2double(tok{1}));
        else
            ord(i) = s(i);
        end
    end
end
end

function trial = compute_trialindex_local(S)
trial = nan(height(S),1);
if ismember('scene_id', S.Properties.VariableNames)
    trial = double(S.scene_id);
end
if ismember('block_id', S.Properties.VariableNames) && ismember('cycle_in_block', S.Properties.VariableNames)
    alt = (double(S.block_id)-1)*6 + double(S.cycle_in_block);
    bad = ~isfinite(trial) | trial<1 | trial>12;
    trial(bad) = alt(bad);
end
bad = ~isfinite(trial) | trial<1 | trial>12;
trial(bad) = nan;
end

function metrics = get_curated_metrics_local(S, cfg)
metrics = {"O_alpha","O_theta","O_beta","F_theta"};
try
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
catch
end
metrics = metrics(ismember(string(metrics), string(S.Properties.VariableNames)));
end

function col = pick_group_col(T, a, b)
col = "";
try
    if ~isempty(T) && istable(T)
        if ismember(a, T.Properties.VariableNames)
            col = string(a);
        elseif ismember(b, T.Properties.VariableNames)
            col = string(b);
        end
    end
catch
    col = "";
end
end

function g = normalize_high_low_local(x)
s = string(x);
s = strtrim(s);
sl = lower(s);
g = repmat("", numel(s), 1);
g(ismember(sl,["high","1","高","h"])) = "High";
g(ismember(sl,["low","0","低","l"])) = "Low";
mask = (s=="High" | s=="Low");
g(mask) = s(mask);
end

function Tsum = summarize_overall_scene(S, cfg)
Tsum = table();
try
    metrics = {"O_alpha","O_theta","O_beta","F_theta"};
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
    rows = {};
    for mi=1:numel(metrics)
        m = string(metrics{mi});
        if ~ismember(m, S.Properties.VariableNames)
            continue;
        end
        y = double(S.(char(m)));
        n = sum(~isnan(y));
        if n==0, continue; end
        mu = mean(y,'omitnan');
        se = std(y,'omitnan')/sqrt(n);
        rows(end+1,:) = {m, n, mu, se}; %#ok<AGROW>
    end
    if ~isempty(rows)
        Tsum = cell2table(rows, 'VariableNames', {'metric','N','mean','sem'});
    end
catch
    Tsum = table();
end
end

function Tsum = summarize_scene_by_experience(S, cfg)
Tsum = table();
try
    metrics = {"O_alpha","O_theta","O_beta","F_theta"};
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
    hasCx = ismember('Complexity', S.Properties.VariableNames);
    rows = {};
    for mi=1:numel(metrics)
        m = string(metrics{mi});
        if ~ismember(m, S.Properties.VariableNames)
            continue;
        end
        y = double(S.(char(m)));
        grp = string(S.ExperienceGroupCurated);
        if hasCx
            cx = string(S.Complexity);
            [G, g, c] = findgroups(grp, cx);
            mu = splitapply(@(x) mean(x,'omitnan'), y, G);
            se = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), y, G);
            n = splitapply(@(x) sum(~isnan(x)), y, G);
            for i=1:numel(mu)
                rows(end+1,:) = {m, g(i), c(i), n(i), mu(i), se(i)}; %#ok<AGROW>
            end
        else
            [G, g] = findgroups(grp);
            mu = splitapply(@(x) mean(x,'omitnan'), y, G);
            se = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), y, G);
            n = splitapply(@(x) sum(~isnan(x)), y, G);
            for i=1:numel(mu)
                rows(end+1,:) = {m, g(i), "", n(i), mu(i), se(i)}; %#ok<AGROW>
            end
        end
    end
    if ~isempty(rows)
        Tsum = cell2table(rows, 'VariableNames', {'metric','experience_group','complexity','N','mean','sem'});
    end
catch
    Tsum = table();
end
end

function Tsum = summarize_overall_by_factors(S, cfg)
Tsum = table();
try
    if ~ismember('WWR', S.Properties.VariableNames) || ~ismember('Complexity', S.Properties.VariableNames)
        return;
    end
    metrics = {"O_alpha","O_theta","O_beta","F_theta"};
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
    rows = {};
    for mi=1:numel(metrics)
        m = string(metrics{mi});
        if ~ismember(m, S.Properties.VariableNames)
            continue;
        end
        y = double(S.(char(m)));
        w = string(S.WWR);
        c = string(S.Complexity);
        [G, ww, cc] = findgroups(w, c);
        mu = splitapply(@(x) mean(x,'omitnan'), y, G);
        se = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), y, G);
        n = splitapply(@(x) sum(~isnan(x)), y, G);
        for i=1:numel(mu)
            rows(end+1,:) = {m, ww(i), cc(i), n(i), mu(i), se(i)}; %#ok<AGROW>
        end
    end
    if ~isempty(rows)
        Tsum = cell2table(rows, 'VariableNames', {'metric','WWR','Complexity','N','mean','sem'});
    end
catch
    Tsum = table();
end
end

function Tsum = summarize_scene_by_experience_factors(S, cfg)
Tsum = table();
try
    if ~ismember('WWR', S.Properties.VariableNames) || ~ismember('Complexity', S.Properties.VariableNames)
        return;
    end
    metrics = {"O_alpha","O_theta","O_beta","F_theta"};
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
    rows = {};
    for mi=1:numel(metrics)
        m = string(metrics{mi});
        if ~ismember(m, S.Properties.VariableNames)
            continue;
        end
        y = double(S.(char(m)));
        g = string(S.ExperienceGroupCurated);
        w = string(S.WWR);
        c = string(S.Complexity);
        [G, gg, ww, cc] = findgroups(g, w, c);
        mu = splitapply(@(x) mean(x,'omitnan'), y, G);
        se = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), y, G);
        n = splitapply(@(x) sum(~isnan(x)), y, G);
        for i=1:numel(mu)
            rows(end+1,:) = {m, gg(i), ww(i), cc(i), n(i), mu(i), se(i)}; %#ok<AGROW>
        end
    end
    if ~isempty(rows)
        Tsum = cell2table(rows, 'VariableNames', {'metric','experience_group','WWR','Complexity','N','mean','sem'});
    end
catch
    Tsum = table();
end
end

function Trec = summarize_recovery_by_experience(P)
Trec = table();
try
    if ~ismember('delta_O_alpha', P.Properties.VariableNames)
        return;
    end
    grp = string(P.ExperienceGroupCurated);
    y = double(P.delta_O_alpha);
    [G, g] = findgroups(grp);
    mu = splitapply(@(x) mean(x,'omitnan'), y, G);
    se = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), y, G);
    n = splitapply(@(x) sum(~isnan(x)), y, G);
    Trec = table(g, n, mu, se, 'VariableNames', {'experience_group','N','mean_delta_O_alpha','sem_delta_O_alpha'});
catch
    Trec = table();
end
end

function copy_matching_figs(srcDir, dstDir, keyword)
if ~exist(srcDir,'dir'); return; end
D = dir(fullfile(srcDir, '*.png'));
for i=1:numel(D)
    if contains(lower(D(i).name), lower(keyword))
        copyfile(fullfile(D(i).folder, D(i).name), fullfile(dstDir, D(i).name), 'f');
    end
end
end

function mirror_task_branch(taskRoot, branchName, dstRoot)
if ~exist(taskRoot,'dir'); return; end
for tag = ["raw","qc"]
    candTbl = {
        fullfile(taskRoot, char(tag), 'tables', branchName), ...
        fullfile(taskRoot, 'tables', char(tag), branchName) ...
    };
    candFig = {
        fullfile(taskRoot, char(tag), 'figures', branchName), ...
        fullfile(taskRoot, 'figures', char(tag), branchName) ...
    };
    candRep = {
        fullfile(taskRoot, char(tag), 'reports', branchName), ...
        fullfile(taskRoot, 'reports', char(tag), branchName) ...
    };

    fpTag = fullfile(dstRoot, char(tag));
    fpTbl = fullfile(fpTag, 'tables');
    fpFig = fullfile(fpTag, 'figures');
    fpRep = fullfile(fpTag, 'report');
    mkdir_if_needed(fpTag); mkdir_if_needed(fpTbl); mkdir_if_needed(fpFig); mkdir_if_needed(fpRep);

    for i=1:numel(candTbl)
        if exist(candTbl{i},'dir')
            copy_dir_contents(candTbl{i}, fpTbl);
            break;
        end
    end
    for i=1:numel(candFig)
        if exist(candFig{i},'dir')
            copy_dir_contents(candFig{i}, fpFig);
            break;
        end
    end
    for i=1:numel(candRep)
        if exist(candRep{i},'dir')
            copy_dir_contents(candRep{i}, fpRep);
            break;
        end
    end
end
end

function copy_dir_contents(srcDir, dstDir)
D = dir(srcDir);
for i=1:numel(D)
    if any(strcmp(D(i).name,{'.','..'})); continue; end
    src = fullfile(D(i).folder, D(i).name);
    dst = fullfile(dstDir, D(i).name);
    if D(i).isdir
        mkdir_if_needed(dst);
        copy_dir_contents(src, dst);
    else
        copyfile(src, dst, 'f');
    end
end
end

function build_inferential_overall_summary(fp_sum, fp_tbl, fp_fig, fp_rep)
files = {
    fullfile(fp_sum,'analysis','MASTER_REPORT.md'), ...
    fullfile(fp_sum,'analysis-2','MASTER_REPORT.md') ...
};
for i=1:numel(files)
    if exist(files{i},'file')
        copyfile(files{i}, fullfile(fp_rep,'MASTER_REPORT.md'), 'f');
        break;
    end
end

% file-level index for overall-oriented task outputs
rows = {};
cands = {
    fullfile(fp_sum,'analysis','task4_core_lmm_suite','raw'), ...
    fullfile(fp_sum,'analysis','task4_core_lmm_suite','qc'), ...
    fullfile(fp_sum,'analysis','task5_peakindex_invertedu','raw'), ...
    fullfile(fp_sum,'analysis','task5_peakindex_invertedu','qc') ...
};
for ci=1:numel(cands)
    d = cands{ci};
    if ~exist(d,'dir'); continue; end
    F = dir(fullfile(d,'**','*.csv'));
    for i=1:numel(F)
        rows(end+1,:) = {string(F(i).name), string(fullfile(F(i).folder, F(i).name))}; %#ok<AGROW>
    end
end
if ~isempty(rows)
    T = cell2table(rows, 'VariableNames', {'file_name','source_path'});
    writetable(T, fullfile(fp_tbl,'overall_inferential_file_index.csv'));
    try
        task = repmat("", height(T), 1);
        for i=1:height(T)
            s = lower(string(T.source_path(i)));
            if contains(s,'task4_core_lmm_suite')
                task(i) = 'task4_core_lmm_suite';
            elseif contains(s,'task5_peakindex_invertedu')
                task(i) = 'task5_peakindex_invertedu';
            else
                task(i) = 'other';
            end
        end
        T.task = task;
        T2 = groupsummary(T, 'task', 'numel', 'file_name');
        if ~isempty(T2)
            writetable(T2, fullfile(fp_tbl,'overall_inferential_task_counts.csv'));
            plot_task_count_bar(T2, fullfile(fp_fig, 'overall_inferential_task_counts.png'));
        end
    catch
    end
end

fid = fopen(fullfile(fp_rep,'README.md'),'w');
if fid>0
    fprintf(fid, '# Inferential / Overall\n\n');
    fprintf(fid, 'This folder is the curated entry point for full-sample inferential outputs.\n\n');
    fprintf(fid, 'Current repository state still stores detailed inferential machinery in task-based folders under `batch/analysis/`.\n');
    fprintf(fid, 'If present, `MASTER_REPORT.md` is copied here as the fastest overall inferential summary.\n');
    fprintf(fid, 'Additional file-level maps may be provided via `overall_inferential_file_index.csv` and `overall_inferential_task_counts.csv`.\n');
    fclose(fid);
end
end

function plot_overall_metric_bar(T, tag, fp_fig, cfg)
if isempty(T) || height(T)==0
    return;
end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 920 430]);
ax = axes(fig); hold(ax,'on'); style_axes(ax);
M = string(T.metric);
vals = double(T.mean);
sem = double(T.sem);
x = 1:numel(vals);
bar(ax, x, vals, 0.58, 'FaceColor',[0.76 0.86 0.94], 'EdgeColor','none');
errorbar(ax, x, vals, sem, 'Color',[0.15 0.15 0.15], 'LineStyle','none', 'LineWidth',1.2, 'CapSize',10);
set(ax,'XTick',x,'XTickLabel',cellstr(M));
ylabel(ax,'Mean ± SEM');
title(ax, sprintf('Overall core metric means [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
for i=1:numel(x)
    text(ax, x(i), vals(i)+max(0.01,sem(i)*1.25), sprintf('%.3f\nN=%d', vals(i), round(T.N(i))), ...
        'HorizontalAlignment','center', 'FontSize',7.5, 'Color',[0.18 0.18 0.18], 'BackgroundColor',[1 1 1], 'Margin',1.5);
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('overall_metric_bar_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function plot_experience_metric_bar(T, tag, fp_fig, cfg)
if isempty(T) || height(T)==0
    return;
end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[80 80 1200 900]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Experience descriptive core metric means [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
metrics = unique(string(T.metric), 'stable');
for mi=1:min(numel(metrics),4)
    m = metrics(mi);
    nexttile; hold on; style_axes(gca);
    X = T(string(T.metric)==m,:);
    groups = unique(string(X.experience_group), 'stable');
    colors = [0.18 0.49 0.72; 0.88 0.47 0.18];
    vals = nan(numel(groups),1); sem = vals; ns = vals;
    for gi=1:numel(groups)
        idx = string(X.experience_group)==groups(gi);
        vals(gi) = double(X.mean(find(idx,1,'first')));
        sem(gi) = double(X.sem(find(idx,1,'first')));
        ns(gi) = double(X.N(find(idx,1,'first')));
    end
    x = 1:numel(groups);
    for gi=1:numel(groups)
        bar(x(gi), vals(gi), 0.58, 'FaceColor', colors(min(gi,2),:), 'EdgeColor','none');
    end
    errorbar(x, vals, sem, 'Color',[0.15 0.15 0.15], 'LineStyle','none', 'LineWidth',1.2, 'CapSize',10);
    set(gca,'XTick',x,'XTickLabel',cellstr(groups));
    title(strrep(char(m),'_','\_'),'Interpreter','none');
    ylabel('Mean ± SEM');
    for gi=1:numel(groups)
        text(x(gi), vals(gi)+max(0.01,sem(gi)*1.25), sprintf('%.3f\nN=%d', vals(gi), round(ns(gi))), ...
            'HorizontalAlignment','center', 'FontSize',7.5, 'Color',[0.18 0.18 0.18], 'BackgroundColor',[1 1 1], 'Margin',1.5);
    end
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('experience_metric_bar_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function Tover = summarize_recovery_overall(P)
Tover = table();
try
    metrics = intersect({'delta_O_alpha','delta_O_theta','delta_O_beta'}, P.Properties.VariableNames, 'stable');
    if isempty(metrics)
        return;
    end
    rows = {};
    for mi=1:numel(metrics)
        m = string(metrics{mi});
        y = double(P.(char(m)));
        n = sum(~isnan(y));
        if n==0, continue; end
        mu = mean(y,'omitnan');
        se = std(y,'omitnan')/sqrt(n);
        rows(end+1,:) = {m, n, mu, se}; %#ok<AGROW>
    end
    if ~isempty(rows)
        Tover = cell2table(rows, 'VariableNames', {'metric','N','mean','sem'});
    end
catch
    Tover = table();
end
end

function plot_recovery_overall(T, tag, fp_fig, cfg)
if isempty(T) || height(T)==0
    return;
end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 900 420]);
ax = axes(fig); hold(ax,'on'); style_axes(ax);
M = string(T.metric);
vals = double(T.mean); sem = double(T.sem);
x = 1:numel(vals);
bar(ax, x, vals, 0.58, 'FaceColor',[0.79 0.89 0.81], 'EdgeColor','none');
errorbar(ax, x, vals, sem, 'Color',[0.15 0.15 0.15], 'LineStyle','none', 'LineWidth',1.2, 'CapSize',10);
yline(ax, 0, '--', 'Color',[0.4 0.4 0.4]);
set(ax,'XTick',x,'XTickLabel',cellstr(M));
ylabel(ax,'Mean ± SEM');
title(ax, sprintf('Overall recovery summary [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
for i=1:numel(x)
    text(ax, x(i), vals(i)+sign(vals(i)+eps)*max(0.01,sem(i)*1.25), sprintf('%.3f\nN=%d', vals(i), round(T.N(i))), ...
        'HorizontalAlignment','center', 'FontSize',7.5, 'Color',[0.18 0.18 0.18], 'BackgroundColor',[1 1 1], 'Margin',1.5);
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('overall_recovery_bar_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function plot_recovery_experience(T, tag, fp_fig, cfg)
if isempty(T) || height(T)==0
    return;
end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 700 420]);
ax = axes(fig); hold(ax,'on'); style_axes(ax);
groups = string(T.experience_group);
vals = double(T.mean_delta_O_alpha); sem = double(T.sem_delta_O_alpha);
x = 1:numel(vals);
colors = [0.18 0.49 0.72; 0.88 0.47 0.18];
for i=1:numel(x)
    bar(x(i), vals(i), 0.58, 'FaceColor', colors(min(i,2),:), 'EdgeColor','none');
end
errorbar(ax, x, vals, sem, 'Color',[0.15 0.15 0.15], 'LineStyle','none', 'LineWidth',1.2, 'CapSize',10);
yline(ax, 0, '--', 'Color',[0.4 0.4 0.4]);
set(ax,'XTick',x,'XTickLabel',cellstr(groups));
ylabel(ax,'Mean delta_O_alpha ± SEM');
title(ax, sprintf('Experience recovery summary [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
for i=1:numel(x)
    text(ax, x(i), vals(i)+sign(vals(i)+eps)*max(0.01,sem(i)*1.25), sprintf('%.3f\nN=%d', vals(i), round(T.N(i))), ...
        'HorizontalAlignment','center', 'FontSize',7.5, 'Color',[0.18 0.18 0.18], 'BackgroundColor',[1 1 1], 'Margin',1.5);
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('experience_recovery_bar_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function plot_overall_factor_grid(T, tag, fp_fig, cfg)
if isempty(T) || height(T)==0
    return;
end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[80 80 1200 900]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Overall descriptive summary by WWR × Complexity [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
metrics = unique(string(T.metric), 'stable');
for mi=1:min(numel(metrics),4)
    m = metrics(mi);
    nexttile;
    hold on;
    style_axes(gca);
    X = T(string(T.metric)==m,:);
    plot_factor_lines(X, false);
    title(strrep(char(m),'_','\_'),'Interpreter','none');
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('overall_factor_grid_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function plot_experience_factor_grid(T, tag, fp_fig, cfg)
if isempty(T) || height(T)==0
    return;
end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[80 80 1200 900]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Experience descriptive summary by WWR [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
metrics = unique(string(T.metric), 'stable');
for mi=1:min(numel(metrics),4)
    m = metrics(mi);
    nexttile;
    hold on;
    style_axes(gca);
    X = T(string(T.metric)==m,:);
    plot_experience_lines_collapsed_complexity(X);
    title(strrep(char(m),'_','\_'),'Interpreter','none');
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('experience_factor_grid_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function plot_factor_lines(X, hasGroup)
colors = [0.18 0.49 0.72; 0.88 0.47 0.18; 0.20 0.62 0.52; 0.62 0.42 0.78];
wwrLevels = [15 45 75];
if hasGroup
    groups = unique(string(X.experience_group), 'stable');
    if isempty(groups), groups = ["Low","High"]; end
    cxLevels = unique(string(X.Complexity), 'stable');
    if isempty(cxLevels), cxLevels = ["0","1"]; end
    ci = 0;
    for gi=1:numel(groups)
        for cxi=1:numel(cxLevels)
            ci = ci + 1;
            g = groups(gi); cx = cxLevels(cxi);
            Ti = X(string(X.experience_group)==g & string(X.Complexity)==cx,:);
            if isempty(Ti), continue; end
            [xpos, ypos, ysem, yn] = map_factor_points_to_discrete_wwr(Ti, wwrLevels);
            ls = '-'; if contains(lower(cx),'1') || contains(lower(cx),'high'), ls='--'; end
            errorbar(xpos, ypos, ysem, 'o', 'LineWidth',1.5, 'Color', colors(min(ci,size(colors,1)),:), ...
                'MarkerFaceColor', colors(min(ci,size(colors,1)),:), 'DisplayName', sprintf('%s | %s', g, cx));
            plot(xpos, ypos, ls, 'LineWidth',1.7, 'Color', colors(min(ci,size(colors,1)),:), 'HandleVisibility','off');
            annotate_points(xpos, ypos, ysem, yn);
        end
    end
    legend('Location','best');
else
    cxLevels = unique(string(X.Complexity), 'stable');
    if isempty(cxLevels), cxLevels = ["0","1"]; end
    for cxi=1:numel(cxLevels)
        cx = cxLevels(cxi);
        Ti = X(string(X.Complexity)==cx,:);
        if isempty(Ti), continue; end
        [xpos, ypos, ysem, yn] = map_factor_points_to_discrete_wwr(Ti, wwrLevels);
        ls = '-'; if contains(lower(cx),'1') || contains(lower(cx),'high'), ls='--'; end
        errorbar(xpos, ypos, ysem, 'o', 'LineWidth',1.5, 'Color', colors(cxi,:), ...
            'MarkerFaceColor', colors(cxi,:), 'DisplayName', sprintf('Complexity %s', cx));
        plot(xpos, ypos, ls, 'LineWidth',1.7, 'Color', colors(cxi,:), 'HandleVisibility','off');
        annotate_points(xpos, ypos, ysem, yn);
    end
    legend('Location','best');
end
set(gca, 'XTick', wwrLevels, 'XTickLabel', compose('%d', wwrLevels));
xlim([12 78]);
xlabel('WWR');
ylabel('Mean ± SEM');
end

function [xpos, ypos, ysem, yn] = map_factor_points_to_discrete_wwr(Ti, wwrLevels)
TiWWR = str2double(string(Ti.WWR));
ypos = nan(size(wwrLevels));
ysem = nan(size(wwrLevels));
yn = nan(size(wwrLevels));
for i = 1:numel(wwrLevels)
    idx = find(TiWWR == wwrLevels(i), 1, 'first');
    if ~isempty(idx)
        ypos(i) = double(Ti.mean(idx));
        ysem(i) = double(Ti.sem(idx));
        yn(i) = double(Ti.N(idx));
    end
end
xpos = wwrLevels;
end

function plot_experience_lines_collapsed_complexity(X)
wwrLevels = [15 45 75];
colors = [0.18 0.49 0.72; 0.88 0.47 0.18];
groups = unique(string(X.experience_group), 'stable');
if isempty(groups), groups = ["Low","High"]; end
for gi = 1:numel(groups)
    g = groups(gi);
    Ti = X(string(X.experience_group)==g,:);
    if isempty(Ti), continue; end
    [xpos, ypos, ysem, yn] = collapse_group_over_complexity(Ti, wwrLevels);
    errorbar(xpos, ypos, ysem, 'o-', 'LineWidth',1.6, 'Color', colors(min(gi,2),:), ...
        'MarkerFaceColor', colors(min(gi,2),:), 'DisplayName', char(g));
    annotate_points(xpos, ypos, ysem, yn);
end
set(gca, 'XTick', wwrLevels, 'XTickLabel', compose('%d', wwrLevels));
xlim([12 78]);
xlabel('WWR');
ylabel('Mean ± SEM');
legend('Location','best');
end

function [xpos, ypos, ysem, yn] = collapse_group_over_complexity(Ti, wwrLevels)
TiWWR = str2double(string(Ti.WWR));
ypos = nan(size(wwrLevels));
ysem = nan(size(wwrLevels));
yn = nan(size(wwrLevels));
for i = 1:numel(wwrLevels)
    idx = find(TiWWR == wwrLevels(i));
    if isempty(idx), continue; end
    m = double(Ti.mean(idx));
    s = double(Ti.sem(idx));
    n = double(Ti.N(idx));
    ypos(i) = mean(m, 'omitnan');
    ysem(i) = mean(s, 'omitnan');
    yn(i) = sum(n, 'omitnan');
end
xpos = wwrLevels;
end

function annotate_points(x, y, se, n)
for i=1:numel(x)
    if ~isfinite(y(i)), continue; end
    bump = max(0.01, abs(se(i))*1.25 + 0.012*max(1,abs(y(i))));
    % stagger labels slightly to reduce overlap in line plots
    if mod(i,2)==0
        yy = y(i) + 1.55*bump;
    else
        yy = y(i) + bump;
    end
    text(x(i), yy, sprintf('%.3f\nN=%d', y(i), round(n(i))), ...
        'HorizontalAlignment','center', 'FontSize',7.5, 'Color',[0.18 0.18 0.18], ...
        'BackgroundColor',[1 1 1], 'Margin',1.5);
end
end

function style_axes(ax)
set(ax, 'Box','off', 'LineWidth',0.8, 'FontName','Times New Roman', 'FontSize',11, 'Color','w');
grid(ax, 'on');
ax.GridAlpha = 0.08;
ax.GridColor = [0 0 0];
ax.XColor = [0.18 0.18 0.18];
ax.YColor = [0.18 0.18 0.18];
end

function dpi = get_dpi(cfg)
dpi = 300;
try
    if isfield(cfg,'paper_dpi') && ~isempty(cfg.paper_dpi)
        dpi = double(cfg.paper_dpi);
    end
catch
end
end

function build_inferential_overall_models(S, cfg, tag, fp_tbl, fp_fig, fp_rep)
if exist('fitlme','file') ~= 2
    return;
end
if ~ismember('subject_id', S.Properties.VariableNames) || ~ismember('WWR', S.Properties.VariableNames) || ~ismember('Complexity', S.Properties.VariableNames)
    return;
end
metrics = {"O_alpha","O_theta","O_beta","F_theta"};
try
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
catch
end
rows = {};
for mi=1:numel(metrics)
    m = string(metrics{mi});
    if ~ismember(m, S.Properties.VariableNames), continue; end
    T = table();
    T.Subject = categorical(string(S.subject_id));
    T.WWR = categorical(normalize_wwr_local(S.WWR), {'15','45','75'});
    T.Complexity = categorical(normalize_complexity_local(S.Complexity), {'ComplexityLow','ComplexityHigh'});
    T.EEG = double(S.(char(m)));
    keep = ~isnan(T.EEG) & ~isundefined(T.WWR) & ~isundefined(T.Complexity);
    T = T(keep,:);
    if height(T) < 20, continue; end
    try
        lme = fitlme(T, 'EEG ~ WWR*Complexity + (1|Subject)');
        A = to_table_compat(anova(lme,'DFMethod','Satterthwaite'));
        writetable(A, fullfile(fp_tbl, sprintf('overall_lmm_anova_%s_%s.csv', lower(char(m)), tag)));
        pWWR = term_p(A, 'WWR');
        pCx = term_p(A, 'Complexity');
        pWxC = interaction_p(A, 'WWR', 'Complexity');
        fWWR = term_f(A, 'WWR');
        fCx = term_f(A, 'Complexity');
        fWxC = interaction_f(A, 'WWR', 'Complexity');
        eWWR = term_eta(A, 'WWR');
        eCx = term_eta(A, 'Complexity');
        eWxC = interaction_eta(A, 'WWR', 'Complexity');
        rows(end+1,:) = {m, height(T), numel(unique(string(T.Subject))), pWWR, pCx, pWxC, fWWR, fCx, fWxC, eWWR, eCx, eWxC}; %#ok<AGROW>
    catch
    end
end
if isempty(rows), return; end
Tout = cell2table(rows, 'VariableNames', {'metric','n_rows','n_subjects','p_WWR','p_Complexity','p_WWRxComplexity','F_WWR','F_Complexity','F_WWRxComplexity','eta_WWR','eta_Complexity','eta_WWRxComplexity'});
writetable(Tout, fullfile(fp_tbl, sprintf('overall_inferential_summary_%s.csv', tag)));
plot_inferential_overall_heatmap(Tout, tag, fp_fig, cfg);
write_inferential_overall_readme(fp_rep, tag);
end

function build_inferential_experience_models(S, cfg, tag, fp_tbl, fp_fig, fp_rep)
if exist('fitlme','file') ~= 2
    return;
end
expCol = pick_group_col(S, 'ExperienceGroup', 'Experience');
if strlength(expCol)==0 || ~ismember('subject_id', S.Properties.VariableNames) || ~ismember('WWR', S.Properties.VariableNames) || ~ismember('Complexity', S.Properties.VariableNames)
    return;
end
metrics = {"O_alpha","O_theta","O_beta","F_theta"};
try
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
catch
end
rows = {};
for mi=1:numel(metrics)
    m = string(metrics{mi});
    if ~ismember(m, S.Properties.VariableNames), continue; end
    T = table();
    T.Subject = categorical(string(S.subject_id));
    T.WWR = categorical(normalize_wwr_local(S.WWR), {'15','45','75'});
    T.Complexity = categorical(normalize_complexity_local(S.Complexity), {'ComplexityLow','ComplexityHigh'});
    T.Group = categorical(normalize_high_low_local(S.(char(expCol))), {'Low','High'});
    T.EEG = double(S.(char(m)));
    keep = ~isnan(T.EEG) & ~isundefined(T.WWR) & ~isundefined(T.Complexity) & ~isundefined(T.Group);
    T = T(keep,:);
    if height(T) < 20, continue; end
    try
        lme = fitlme(T, 'EEG ~ WWR*Complexity*Group + (1|Subject)');
        A = to_table_compat(anova(lme,'DFMethod','Satterthwaite'));
        writetable(A, fullfile(fp_tbl, sprintf('experience_lmm_anova_%s_%s.csv', lower(char(m)), tag)));
        rows(end+1,:) = {m, height(T), numel(unique(string(T.Subject))), ...
            term_p(A,'WWR'), term_p(A,'Complexity'), term_p(A,'Group'), ...
            interaction_p(A,'WWR','Complexity'), interaction_p(A,'WWR','Group'), interaction_p(A,'Complexity','Group'), ...
            threeway_p(A,'WWR','Complexity','Group'), ...
            term_f(A,'WWR'), term_f(A,'Complexity'), term_f(A,'Group'), ...
            interaction_f(A,'WWR','Complexity'), interaction_f(A,'WWR','Group'), interaction_f(A,'Complexity','Group'), ...
            threeway_f(A,'WWR','Complexity','Group'), ...
            term_eta(A,'WWR'), term_eta(A,'Complexity'), term_eta(A,'Group'), ...
            interaction_eta(A,'WWR','Complexity'), interaction_eta(A,'WWR','Group'), interaction_eta(A,'Complexity','Group'), ...
            threeway_eta(A,'WWR','Complexity','Group')}; %#ok<AGROW>
    catch
    end
end
if isempty(rows), return; end
Tout = cell2table(rows, 'VariableNames', {'metric','n_rows','n_subjects','p_WWR','p_Complexity','p_Group','p_WWRxComplexity','p_WWRxGroup','p_ComplexityxGroup','p_threeway','F_WWR','F_Complexity','F_Group','F_WWRxComplexity','F_WWRxGroup','F_ComplexityxGroup','F_threeway','eta_WWR','eta_Complexity','eta_Group','eta_WWRxComplexity','eta_WWRxGroup','eta_ComplexityxGroup','eta_threeway'});
writetable(Tout, fullfile(fp_tbl, sprintf('experience_inferential_summary_%s.csv', tag)));
plot_inferential_experience_heatmap(Tout, tag, fp_fig, cfg);
write_inferential_experience_readme(fp_rep, tag);
end

function build_inferential_overall_recovery(P, tag, fp_tbl, fp_fig, fp_rep)
if isempty(P) || ~istable(P)
    return;
end
metrics = intersect({'delta_O_alpha','delta_O_theta','delta_O_beta'}, P.Properties.VariableNames, 'stable');
if isempty(metrics)
    return;
end
rows = {};
for mi=1:numel(metrics)
    m = string(metrics{mi});
    y = double(P.(char(m)));
    y = y(isfinite(y));
    n = numel(y);
    if n < 3, continue; end
    mu = mean(y,'omitnan');
    se = std(y,'omitnan')/sqrt(n);
    p = NaN; tstat = NaN; df = NaN; dz = NaN;
    try
        [~, p, ~, st] = ttest(y, 0);
        tstat = st.tstat; df = st.df;
        dz = mu / std(y,'omitnan');
    catch
    end
    rows(end+1,:) = {m, n, mu, se, p, tstat, df, dz}; %#ok<AGROW>
end
if isempty(rows), return; end
T = cell2table(rows, 'VariableNames', {'metric','N','mean','sem','p_value','t_stat','df','cohen_dz'});
writetable(T, fullfile(fp_tbl, sprintf('overall_recovery_inferential_summary_%s.csv', tag)));
plot_recovery_inferential_overall(T, tag, fp_fig);
end

function build_inferential_experience_recovery(P, tag, fp_tbl, fp_fig, fp_rep)
if isempty(P) || ~istable(P)
    return;
end
expCol = pick_group_col(P, 'ExperienceGroup', 'Experience');
if strlength(expCol)==0
    return;
end
metrics = intersect({'delta_O_alpha','delta_O_theta','delta_O_beta'}, P.Properties.VariableNames, 'stable');
if isempty(metrics)
    return;
end
rows = {};
grp = normalize_high_low_local(P.(char(expCol)));
for mi=1:numel(metrics)
    m = string(metrics{mi});
    y = double(P.(char(m)));
    keep = isfinite(y) & strlength(grp)>0;
    y = y(keep); g = grp(keep);
    if numel(y) < 6, continue; end
    yL = y(g=="Low"); yH = y(g=="High");
    if numel(yL) < 3 || numel(yH) < 3, continue; end
    p = NaN; tstat = NaN; df = NaN; d = NaN;
    try
        [~, p, ~, st] = ttest2(yH, yL);
        tstat = st.tstat; df = st.df;
        sP = sqrt(((numel(yH)-1)*var(yH) + (numel(yL)-1)*var(yL)) / max(1,(numel(yH)+numel(yL)-2)));
        d = (mean(yH,'omitnan') - mean(yL,'omitnan')) / sP;
    catch
    end
    rows(end+1,:) = {m, numel(yH), numel(yL), mean(yH,'omitnan'), mean(yL,'omitnan'), p, tstat, df, d}; %#ok<AGROW>
end
if isempty(rows), return; end
T = cell2table(rows, 'VariableNames', {'metric','N_High','N_Low','mean_High','mean_Low','p_group','t_group','df','cohen_d'});
writetable(T, fullfile(fp_tbl, sprintf('experience_recovery_inferential_summary_%s.csv', tag)));
plot_recovery_inferential_experience(T, tag, fp_fig);
end

function build_inferential_overall_trend(S, cfg, tag, fp_tbl, fp_fig, fp_rep)
if exist('fitlme','file') ~= 2
    return;
end
if ~ismember('subject_id', S.Properties.VariableNames) || ~ismember('WWR', S.Properties.VariableNames) || ~ismember('Complexity', S.Properties.VariableNames)
    return;
end
metrics = {"O_alpha","O_theta","O_beta","F_theta"};
try
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
catch
end
rows = {};
for mi=1:numel(metrics)
    m = string(metrics{mi});
    if ~ismember(m, S.Properties.VariableNames), continue; end
    T = table();
    T.Subject = categorical(string(S.subject_id));
    T.Complexity = categorical(normalize_complexity_local(S.Complexity), {'ComplexityLow','ComplexityHigh'});
    w = str2double(string(normalize_wwr_local(S.WWR)));
    T.WWRc = w - 45;
    T.WWRc2 = T.WWRc.^2;
    T.EEG = double(S.(char(m)));
    keep = ~isnan(T.EEG) & ~isundefined(T.Complexity) & isfinite(T.WWRc) & isfinite(T.WWRc2);
    T = T(keep,:);
    if height(T) < 20, continue; end
    try
        lmeL = fitlme(T, 'EEG ~ WWRc + Complexity + (1|Subject)');
        lmeQ = fitlme(T, 'EEG ~ WWRc + WWRc2 + Complexity + (1|Subject)');
        [bL,pL,bQ,pQ] = extract_trend_stats_local(lmeL, lmeQ);
        verdict = trend_verdict(bL,pL,bQ,pQ);
        rows(end+1,:) = {m, height(T), numel(unique(string(T.Subject))), bL, pL, bQ, pQ, verdict}; %#ok<AGROW>
    catch
    end
end
if isempty(rows), return; end
Tout = cell2table(rows, 'VariableNames', {'metric','n_rows','n_subjects','beta_linear','p_linear','beta_quadratic','p_quadratic','trend_verdict'});
writetable(Tout, fullfile(fp_tbl, sprintf('overall_wwr_trend_summary_%s.csv', tag)));
plot_trend_heatmap(Tout, sprintf('Overall WWR trend [%s]', tag), fullfile(fp_fig, sprintf('overall_wwr_trend_heatmap_%s.png', tag)), cfg);
end

function build_inferential_experience_trend(S, cfg, tag, fp_tbl, fp_fig, fp_rep)
if exist('fitlme','file') ~= 2
    return;
end
expCol = pick_group_col(S, 'ExperienceGroup', 'Experience');
if strlength(expCol)==0 || ~ismember('subject_id', S.Properties.VariableNames) || ~ismember('WWR', S.Properties.VariableNames) || ~ismember('Complexity', S.Properties.VariableNames)
    return;
end
metrics = {"O_alpha","O_theta","O_beta","F_theta"};
try
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
catch
end
rows = {};
for mi=1:numel(metrics)
    m = string(metrics{mi});
    if ~ismember(m, S.Properties.VariableNames), continue; end
    T = table();
    T.Subject = categorical(string(S.subject_id));
    T.Group = categorical(normalize_high_low_local(S.(char(expCol))), {'Low','High'});
    T.Complexity = categorical(normalize_complexity_local(S.Complexity), {'ComplexityLow','ComplexityHigh'});
    w = str2double(string(normalize_wwr_local(S.WWR)));
    T.WWRc = w - 45;
    T.WWRc2 = T.WWRc.^2;
    T.EEG = double(S.(char(m)));
    keep = ~isnan(T.EEG) & ~isundefined(T.Complexity) & ~isundefined(T.Group) & isfinite(T.WWRc) & isfinite(T.WWRc2);
    T = T(keep,:);
    if height(T) < 20, continue; end
    try
        lmeL = fitlme(T, 'EEG ~ WWRc + Complexity + Group + (1|Subject)');
        lmeQ = fitlme(T, 'EEG ~ WWRc + WWRc2 + Complexity + Group + (1|Subject)');
        [bL,pL,bQ,pQ] = extract_trend_stats_local(lmeL, lmeQ);
        verdict = trend_verdict(bL,pL,bQ,pQ);
        rows(end+1,:) = {m, height(T), numel(unique(string(T.Subject))), bL, pL, bQ, pQ, verdict}; %#ok<AGROW>
    catch
    end
end
if isempty(rows), return; end
Tout = cell2table(rows, 'VariableNames', {'metric','n_rows','n_subjects','beta_linear','p_linear','beta_quadratic','p_quadratic','trend_verdict'});
writetable(Tout, fullfile(fp_tbl, sprintf('experience_wwr_trend_summary_%s.csv', tag)));
plot_trend_heatmap(Tout, sprintf('Experience WWR trend [%s]', tag), fullfile(fp_fig, sprintf('experience_wwr_trend_heatmap_%s.png', tag)), cfg);
end

function [bL,pL,bQ,pQ] = extract_trend_stats_local(lmeL, lmeQ)
bL=NaN; pL=NaN; bQ=NaN; pQ=NaN;
try
    C = to_table_compat(lmeL.Coefficients);
    n = string(C.Name); i = find(n=="WWRc",1,'first');
    if ~isempty(i), bL = double(C.Estimate(i)); pL = double(C.pValue(i)); end
    A = to_table_compat(anova(lmeL,'DFMethod','Satterthwaite'));
    pA = term_p(A, 'WWRc'); if isfinite(pA), pL = pA; end
catch
end
try
    C = to_table_compat(lmeQ.Coefficients);
    n = string(C.Name); i = find(n=="WWRc2",1,'first');
    if ~isempty(i), bQ = double(C.Estimate(i)); pQ = double(C.pValue(i)); end
    A = to_table_compat(anova(lmeQ,'DFMethod','Satterthwaite'));
    pA = term_p(A, 'WWRc2'); if isfinite(pA), pQ = pA; end
catch
end
end

function v = trend_verdict(bL,pL,bQ,pQ)
v = "none";
if isfinite(pQ) && pQ < 0.05
    if bQ < 0
        v = "peak_45";
    elseif bQ > 0
        v = "trough_45";
    else
        v = "quadratic";
    end
elseif isfinite(pL) && pL < 0.05
    if bL > 0
        v = "linear_increase";
    elseif bL < 0
        v = "linear_decrease";
    else
        v = "linear";
    end
end
end

function plot_trend_heatmap(T, ttl, fp_png, cfg)
if isempty(T) || height(T)==0, return; end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 980 360]);
ax = axes(fig); hold(ax,'on');
Z = [double(T.p_linear)'; double(T.p_quadratic)'];
V = nan(size(Z)); V(isfinite(Z))=0; V(isfinite(Z) & Z<0.05)=1; V(isfinite(Z) & Z<0.01)=2; V(isfinite(Z) & Z<0.001)=3;
imagesc(ax, V); set(ax,'YDir','normal'); axis(ax,'tight');
colormap(ax, interp1([0 0.5 1], [0.80 0.88 0.94; 0.98 0.98 0.98; 0.22 0.45 0.67], linspace(0,1,256)));
caxis(ax,[0 3]);
set(ax,'XTick',1:height(T),'XTickLabel',cellstr(string(T.metric)),'YTick',1:2,'YTickLabel',{'Linear','Quadratic'});
xtickangle(ax,18); title(ax, ttl, 'Interpreter','none', 'FontWeight','normal');
style_axes(ax);
for c=1:height(T)
    text(ax,c,1,sprintf('p=%.3g%s\nβ=%.3g',double(T.p_linear(c)),star_from_p(double(T.p_linear(c))),double(T.beta_linear(c))), 'HorizontalAlignment','center', 'FontSize',8, 'BackgroundColor',[1 1 1], 'Margin',1);
    text(ax,c,2,sprintf('p=%.3g%s\nβ2=%.3g\n%s',double(T.p_quadratic(c)),star_from_p(double(T.p_quadratic(c))),double(T.beta_quadratic(c)),char(string(T.trend_verdict(c)))), 'HorizontalAlignment','center', 'FontSize',8, 'BackgroundColor',[1 1 1], 'Margin',1);
end
pipeline.export_figure_png(fig, fp_png, get_dpi(cfg));
try; close(fig); catch; end
end

function f = term_f(A, name)
f = NaN;
try
    if all(ismember({'Term'}, A.Properties.VariableNames))
        tt = string(A.Term);
        idx = find(tt==string(name),1,'first');
        if ~isempty(idx), f = get_fstat(A, idx); end
    end
catch
end
end

function f = interaction_f(A, a, b)
f = NaN;
try
    if all(ismember({'Term'}, A.Properties.VariableNames))
        tt = string(A.Term);
        idx = find(contains(tt,a) & contains(tt,b) & contains(tt,':'),1,'first');
        if ~isempty(idx), f = get_fstat(A, idx); end
    end
catch
end
end

function f = threeway_f(A, a, b, c)
f = NaN;
try
    if all(ismember({'Term'}, A.Properties.VariableNames))
        tt = string(A.Term);
        idx = find(contains(tt,a) & contains(tt,b) & contains(tt,c) & contains(tt,':'),1,'first');
        if ~isempty(idx), f = get_fstat(A, idx); end
    end
catch
end
end

function e = term_eta(A, name)
e = NaN;
try
    if all(ismember({'Term'}, A.Properties.VariableNames))
        tt = string(A.Term);
        idx = find(tt==string(name),1,'first');
        if ~isempty(idx), e = get_partial_eta(A, idx); end
    end
catch
end
end

function e = interaction_eta(A, a, b)
e = NaN;
try
    if all(ismember({'Term'}, A.Properties.VariableNames))
        tt = string(A.Term);
        idx = find(contains(tt,a) & contains(tt,b) & contains(tt,':'),1,'first');
        if ~isempty(idx), e = get_partial_eta(A, idx); end
    end
catch
end
end

function e = threeway_eta(A, a, b, c)
e = NaN;
try
    if all(ismember({'Term'}, A.Properties.VariableNames))
        tt = string(A.Term);
        idx = find(contains(tt,a) & contains(tt,b) & contains(tt,c) & contains(tt,':'),1,'first');
        if ~isempty(idx), e = get_partial_eta(A, idx); end
    end
catch
end
end

function f = get_fstat(A, idx)
f = NaN;
vars = string(A.Properties.VariableNames);
for nm = ["FStat","F","FValue","F_stat","Fstat"]
    if any(vars==nm)
        f = double(A.(char(nm))(idx));
        return;
    end
end
end

function e = get_partial_eta(A, idx)
e = NaN;
try
    f = get_fstat(A, idx);
    if ~isfinite(f)
        return;
    end
    vars = string(A.Properties.VariableNames);
    df1 = NaN; df2 = NaN;
    for nm = ["DF1","NumDF","df1"]
        if any(vars==nm)
            df1 = double(A.(char(nm))(idx));
            break;
        end
    end
    for nm = ["DF2","DenDF","df2"]
        if any(vars==nm)
            df2 = double(A.(char(nm))(idx));
            break;
        end
    end
    if isfinite(df1) && isfinite(df2)
        e = (f*df1) / (f*df1 + df2);
    end
catch
end
end

function p = term_p(A, name)
p = NaN;
try
    if all(ismember({'Term','pValue'}, A.Properties.VariableNames))
        tt = string(A.Term);
        idx = find(tt==string(name),1,'first');
        if ~isempty(idx), p = double(A.pValue(idx)); end
    end
catch
end
end

function p = interaction_p(A, a, b)
p = NaN;
try
    if all(ismember({'Term','pValue'}, A.Properties.VariableNames))
        tt = string(A.Term);
        idx = find(contains(tt,a) & contains(tt,b) & contains(tt,':'),1,'first');
        if ~isempty(idx), p = double(A.pValue(idx)); end
    end
catch
end
end

function p = threeway_p(A, a, b, c)
p = NaN;
try
    if all(ismember({'Term','pValue'}, A.Properties.VariableNames))
        tt = string(A.Term);
        idx = find(contains(tt,a) & contains(tt,b) & contains(tt,c) & contains(tt,':'),1,'first');
        if ~isempty(idx), p = double(A.pValue(idx)); end
    end
catch
end
end

function plot_inferential_overall_heatmap(T, tag, fp_fig, cfg)
if isempty(T) || height(T)==0, return; end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 1040 380]);
ax = axes(fig); hold(ax,'on');
Z = [double(T.p_WWR)'; double(T.p_Complexity)'; double(T.p_WWRxComplexity)'];
V = nan(size(Z)); V(isfinite(Z)) = 0; V(isfinite(Z) & Z<0.05)=1; V(isfinite(Z) & Z<0.01)=2; V(isfinite(Z) & Z<0.001)=3;
imagesc(ax, V); set(ax,'YDir','normal'); axis(ax,'tight');
colormap(ax, interp1([0 0.5 1], [0.80 0.88 0.94; 0.98 0.98 0.98; 0.22 0.45 0.67], linspace(0,1,256)));
caxis(ax,[0 3]);
set(ax,'XTick',1:height(T),'XTickLabel',cellstr(string(T.metric)),'YTick',1:3,'YTickLabel',{'WWR','Complexity','WWR×Complexity'});
xtickangle(ax,18); title(ax, sprintf('Inferential / Overall [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
style_axes(ax);
for r=1:size(Z,1)
 for c=1:size(Z,2)
  if isnan(Z(r,c)), continue; end
  text(ax,c,r,sprintf('p=%.3g%s',Z(r,c), star_from_p(Z(r,c))),'HorizontalAlignment','center','FontSize',8, 'BackgroundColor',[1 1 1], 'Margin',1);
 end
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('overall_inferential_heatmap_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function plot_inferential_experience_heatmap(T, tag, fp_fig, cfg)
if isempty(T) || height(T)==0, return; end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 1260 560]);
ax = axes(fig); hold(ax,'on');
Z = [double(T.p_WWR)'; double(T.p_Complexity)'; double(T.p_Group)'; double(T.p_WWRxComplexity)'; double(T.p_WWRxGroup)'; double(T.p_ComplexityxGroup)'; double(T.p_threeway)'];
V = nan(size(Z)); V(isfinite(Z)) = 0; V(isfinite(Z) & Z<0.05)=1; V(isfinite(Z) & Z<0.01)=2; V(isfinite(Z) & Z<0.001)=3;
imagesc(ax, V); set(ax,'YDir','normal'); axis(ax,'tight');
colormap(ax, interp1([0 0.5 1], [0.80 0.88 0.94; 0.98 0.98 0.98; 0.22 0.45 0.67], linspace(0,1,256)));
caxis(ax,[0 3]);
set(ax,'XTick',1:height(T),'XTickLabel',cellstr(string(T.metric)),'YTick',1:7,'YTickLabel',{'WWR','Complexity','Group','WWR×Complexity','WWR×Group','Complexity×Group','3-way'});
xtickangle(ax,18); title(ax, sprintf('Inferential / Experience [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
style_axes(ax);
for r=1:size(Z,1)
 for c=1:size(Z,2)
  if isnan(Z(r,c)), continue; end
  text(ax,c,r,sprintf('p=%.3g%s',Z(r,c), star_from_p(Z(r,c))),'HorizontalAlignment','center','FontSize',8, 'BackgroundColor',[1 1 1], 'Margin',1);
 end
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('experience_inferential_heatmap_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function s = star_from_p(p)
s = '';
if ~isfinite(p), return; end
if p < 0.001
    s = '***';
elseif p < 0.01
    s = '**';
elseif p < 0.05
    s = '*';
end
end

function plot_recovery_inferential_overall(T, tag, fp_fig)
if isempty(T) || height(T)==0, return; end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 920 420]);
ax = axes(fig); hold(ax,'on'); style_axes(ax);
vals = double(T.mean); sem = double(T.sem); x = 1:height(T);
bar(ax, x, vals, 0.58, 'FaceColor',[0.78 0.88 0.82], 'EdgeColor','none');
errorbar(ax, x, vals, sem, 'Color',[0.15 0.15 0.15], 'LineStyle','none', 'LineWidth',1.2, 'CapSize',10);
yline(ax, 0, '--', 'Color',[0.4 0.4 0.4]);
set(ax,'XTick',x,'XTickLabel',cellstr(string(T.metric)));
ylabel(ax,'Mean ± SEM');
title(ax, sprintf('Inferential recovery / Overall [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
for i=1:height(T)
    text(ax, x(i), vals(i)+sign(vals(i)+eps)*max(0.01,sem(i)*1.25), sprintf('p=%.3g%s\nd=%.2f', double(T.p_value(i)), star_from_p(double(T.p_value(i))), double(T.cohen_dz(i))), 'HorizontalAlignment','center', 'FontSize',7.5, 'Color',[0.18 0.18 0.18], 'BackgroundColor',[1 1 1], 'Margin',1.5);
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('overall_recovery_inferential_%s.png', tag)), 300);
try; close(fig); catch; end
end

function plot_recovery_inferential_experience(T, tag, fp_fig)
if isempty(T) || height(T)==0, return; end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 980 420]);
ax = axes(fig); hold(ax,'on'); style_axes(ax);
vals = [double(T.mean_Low), double(T.mean_High)];
sem = zeros(size(vals));
x0 = 1:height(T);
for i=1:height(T)
    bar(ax, x0(i)-0.16, vals(i,1), 0.28, 'FaceColor',[0.18 0.49 0.72], 'EdgeColor','none');
    bar(ax, x0(i)+0.16, vals(i,2), 0.28, 'FaceColor',[0.88 0.47 0.18], 'EdgeColor','none');
    text(ax, x0(i), max(vals(i,:)) + 0.03*max(1,max(vals(:))), sprintf('p=%.3g%s\nd=%.2f', double(T.p_group(i)), star_from_p(double(T.p_group(i))), double(T.cohen_d(i))), 'HorizontalAlignment','center', 'FontSize',7.5, 'BackgroundColor',[1 1 1], 'Margin',1.5);
end
yline(ax, 0, '--', 'Color',[0.4 0.4 0.4]);
set(ax,'XTick',x0,'XTickLabel',cellstr(string(T.metric)));
ylabel(ax,'Group means');
title(ax, sprintf('Inferential recovery / Experience [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
legend({'Low','High'}, 'Location','best');
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('experience_recovery_inferential_%s.png', tag)), 300);
try; close(fig); catch; end
end

function write_inferential_overall_readme(fp_rep, tag)
fid = fopen(fullfile(fp_rep, sprintf('README_%s.md', tag)),'w');
if fid>0
 fprintf(fid,'# Inferential / Overall [%s]\n\n', tag);
 fprintf(fid,'Key summary files:\n');
 fprintf(fid,'- `overall_inferential_summary_%s.csv` (includes p-values, F-statistics, and effect-size proxy columns when available)\n', tag);
 fprintf(fid,'- `overall_inferential_heatmap_%s.png`\n', tag);
 fprintf(fid,'- `overall_wwr_trend_summary_%s.csv`\n', tag);
 fprintf(fid,'- `overall_wwr_trend_heatmap_%s.png`\n', tag);
 fprintf(fid,'- `overall_recovery_inferential_summary_%s.csv`\n', tag);
 fprintf(fid,'- `overall_recovery_inferential_%s.png`\n', tag);
 fprintf(fid,'- `overall_inferential_file_index.csv`\n');
 fprintf(fid,'- `overall_inferential_task_counts.csv`\n');
 fprintf(fid,'- `overall_inferential_task_counts.png`\n');
 fclose(fid);
end
end

function write_inferential_experience_readme(fp_rep, tag)
fid = fopen(fullfile(fp_rep, sprintf('README_%s.md', tag)),'w');
if fid>0
 fprintf(fid,'# Inferential / Experience [%s]\n\n', tag);
 fprintf(fid,'Key summary files:\n');
 fprintf(fid,'- `experience_inferential_summary_%s.csv` (includes p-values, F-statistics, and effect-size proxy columns when available)\n', tag);
 fprintf(fid,'- `experience_inferential_heatmap_%s.png`\n', tag);
 fprintf(fid,'- `experience_wwr_trend_summary_%s.csv`\n', tag);
 fprintf(fid,'- `experience_wwr_trend_heatmap_%s.png`\n', tag);
 fprintf(fid,'- `experience_recovery_inferential_summary_%s.csv`\n', tag);
 fprintf(fid,'- `experience_recovery_inferential_%s.png`\n', tag);
 fprintf(fid,'- mirrored task3 TrialIndex outputs when available (for trial-by-trial neural response significance)\n');
 fprintf(fid,'- `experience_inferential_file_index.csv`\n');
 fprintf(fid,'- `experience_inferential_task_counts.csv`\n');
 fprintf(fid,'- `experience_inferential_task_counts.png`\n');
 fclose(fid);
end
end

function plot_task_count_bar(T, fp_png)
if isempty(T) || height(T)==0
    return;
end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 920 420]);
ax = axes(fig); hold(ax,'on'); style_axes(ax);
vals = double(T.GroupCount);
labels = string(T.task);
x = 1:numel(vals);
bar(ax, x, vals, 0.58, 'FaceColor',[0.82 0.88 0.95], 'EdgeColor','none');
set(ax,'XTick',x,'XTickLabel',cellstr(labels));
xtickangle(ax,20);
ylabel(ax,'Number of CSV files');
title(ax,'Experience inferential task coverage', 'Interpreter','none', 'FontWeight','normal');
for i=1:numel(x)
    text(ax, x(i), vals(i)+max(0.2,0.05*max(vals)), sprintf('%d', round(vals(i))), 'HorizontalAlignment','center', 'FontSize',8, 'Color',[0.18 0.18 0.18], 'BackgroundColor',[1 1 1], 'Margin',1.5);
end
pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function v = to_table_compat(X)
if istable(X), v = X; return; end
try
 if isa(X,'dataset'), v = dataset2table(X); return; end
catch
end
try
 v = struct2table(X); return;
catch
end
error('Unsupported output type');
end

function w = normalize_wwr_local(x)
w = string(x); w = strtrim(w);
for i=1:numel(w)
 tok = regexp(char(w(i)), '(\d+)', 'tokens', 'once');
 if ~isempty(tok), w(i) = string(str2double(tok{1})); end
end
ok = ismember(w,["15","45","75"]); w(~ok) = "";
end

function c = normalize_complexity_local(x)
c = string(x); c = strtrim(c); cl = lower(c); out = repmat("", numel(c), 1);
out(ismember(cl,["low","0","c0","complexitylow"])) = "ComplexityLow";
out(ismember(cl,["high","1","c1","complexityhigh"])) = "ComplexityHigh";
out(c=="ComplexityLow") = "ComplexityLow"; out(c=="ComplexityHigh") = "ComplexityHigh";
isNum = ~isnan(str2double(c)); out(isNum & str2double(c)==0) = "ComplexityLow"; out(isNum & str2double(c)==1) = "ComplexityHigh";
c = out;
end

function build_inferential_experience_summary(fp_sum, fp_tbl, fp_fig, fp_rep)
rows = {};
% Task4 model1/model2/model3 and Task3 key ANOVA summaries for experience branch
cands = {
    fullfile(fp_sum,'analysis','task4_core_lmm_suite','raw','tables','factor_WWR','experience'), ...
    fullfile(fp_sum,'analysis','task4_core_lmm_suite','qc','tables','factor_WWR','experience'), ...
    fullfile(fp_sum,'analysis','task3_trialindex_lmm','raw','tables','experience'), ...
    fullfile(fp_sum,'analysis','task3_trialindex_lmm','qc','tables','experience'), ...
    fullfile(fp_sum,'analysis','task5_peakindex_invertedu','raw','tables','experience'), ...
    fullfile(fp_sum,'analysis','task5_peakindex_invertedu','qc','tables','experience'), ...
    fullfile(fp_sum,'analysis','task6_coremetric_special','raw','tables','experience'), ...
    fullfile(fp_sum,'analysis','task6_coremetric_special','qc','tables','experience'), ...
    fullfile(fp_sum,'analysis','task7_individual_checks','raw','tables','experience'), ...
    fullfile(fp_sum,'analysis','task7_individual_checks','qc','tables','experience') ...
};
for ci=1:numel(cands)
    d = cands{ci};
    if ~exist(d,'dir'); continue; end
    F = dir(fullfile(d,'*.csv'));
    for i=1:numel(F)
        rows(end+1,:) = {string(F(i).name), string(fullfile(F(i).folder, F(i).name))}; %#ok<AGROW>
    end
end
if ~isempty(rows)
    T = cell2table(rows, 'VariableNames', {'file_name','source_path'});
    writetable(T, fullfile(fp_tbl,'experience_inferential_file_index.csv'));

    try
        task = repmat("", height(T), 1);
        for i=1:height(T)
            s = lower(string(T.source_path(i)));
            if contains(s,'task3_trialindex_lmm')
                task(i) = 'task3_trialindex_lmm';
            elseif contains(s,'task4_core_lmm_suite')
                task(i) = 'task4_core_lmm_suite';
            elseif contains(s,'task5_peakindex_invertedu')
                task(i) = 'task5_peakindex_invertedu';
            elseif contains(s,'task6_coremetric_special')
                task(i) = 'task6_coremetric_special';
            elseif contains(s,'task7_individual_checks')
                task(i) = 'task7_individual_checks';
            else
                task(i) = 'other';
            end
        end
        T.task = task;
        T2 = groupsummary(T, 'task', 'numel', 'file_name');
        if ~isempty(T2)
            writetable(T2, fullfile(fp_tbl, 'experience_inferential_task_counts.csv'));
            plot_task_count_bar(T2, fullfile(fp_fig, 'experience_inferential_task_counts.png'));
        end
    catch
    end
end
fid = fopen(fullfile(fp_rep,'README.md'),'w');
if fid>0
    fprintf(fid, '# Inferential / Experience\n\n');
    fprintf(fid, 'This folder is the curated entry point for experience-group inferential outputs.\n');
    fprintf(fid, 'Use `experience_inferential_file_index.csv` as the fastest file-level map into the detailed analysis outputs.\n');
    fprintf(fid, 'Use `experience_inferential_task_counts.csv` as a quick summary of which task families currently contribute files here.\n\n');
    fprintf(fid, 'Primary source tasks currently include:\n');
    fprintf(fid, '- task3_trialindex_lmm\n');
    fprintf(fid, '- task4_core_lmm_suite\n');
    fprintf(fid, '- task5_peakindex_invertedu\n');
    fprintf(fid, '- task6_coremetric_special\n');
    fprintf(fid, '- task7_individual_checks\n\n');
    fprintf(fid, 'SportFreq outputs remain available in the detailed analysis tree, but are not elevated to a first-class curated branch in the redesigned main output surface.\n');
    fclose(fid);
end
end
