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

% Root-level guide for the redesigned main output surface
try
    fid = fopen(fullfile(fp_sum, 'README_CURATED_MAIN.md'),'w');
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
        fprintf(fid, '- `pairs_overall_raw.csv`: full-sample recovery/pair descriptive table\n');
        fprintf(fid, '- `pairs_overall_qc.csv`: QC-filtered recovery/pair descriptive table\n');
        fprintf(fid, '- `overall_recovery_means_*.csv`: overall recovery summary table\n');
        fprintf(fid, '- `overall_metric_bar_*.png`: overall core metric bar figure\n');
        fprintf(fid, '- `overall_recovery_bar_*.png`: overall recovery bar figure\n');
        fprintf(fid, '- `overall_factor_grid_*.png`: overall WWR × Complexity descriptive grid\n');
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
        fprintf(fid, '- `pairs_experience_*.csv`: recovery tables with valid experience labels\n');
        fprintf(fid, '- `experience_recovery_means_*.csv`: recovery summaries by Experience group\n');
        fprintf(fid, '- `experience_metric_bar_*.png`: experience-group core metric bar figure\n');
        fprintf(fid, '- `experience_recovery_bar_*.png`: experience-group recovery bar figure\n');
        fprintf(fid, '- `experience_factor_grid_*.png`: experience-group WWR × Complexity descriptive grid\n');
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

end

function mkdir_if_needed(p)
if ~exist(p,'dir'); mkdir(p); end
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
fid = fopen(fullfile(fp_rep,'README.md'),'w');
if fid>0
    fprintf(fid, '# Inferential / Overall\n\n');
    fprintf(fid, 'This folder is the curated entry point for full-sample inferential outputs.\n\n');
    fprintf(fid, 'Current repository state still stores detailed inferential machinery in task-based folders under `batch/analysis/`.\n');
    fprintf(fid, 'If present, `MASTER_REPORT.md` is copied here as the fastest overall inferential summary.\n');
    fclose(fid);
end
end

function plot_overall_metric_bar(T, tag, fp_fig, cfg)
if isempty(T) || height(T)==0
    return;
end
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 900 420]);
ax = axes(fig); hold(ax,'on'); style_axes(ax);
M = string(T.metric);
vals = double(T.mean);
sem = double(T.sem);
x = 1:numel(vals);
b = bar(ax, x, vals, 0.62, 'FaceColor',[0.78 0.86 0.94], 'EdgeColor','none'); %#ok<NASGU>
errorbar(ax, x, vals, sem, 'k.', 'LineWidth',1.4);
set(ax,'XTick',x,'XTickLabel',cellstr(M));
ylabel(ax,'Mean ± SEM');
title(ax, sprintf('Overall core metric means [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
for i=1:numel(x)
    text(ax, x(i), vals(i)+max(0.01,sem(i)*1.2), sprintf('%.3f\nN=%d', vals(i), round(T.N(i))), 'HorizontalAlignment','center', 'FontSize',8, 'Color',[0.2 0.2 0.2]);
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
        bar(x(gi), vals(gi), 0.62, 'FaceColor', colors(min(gi,2),:), 'EdgeColor','none');
    end
    errorbar(x, vals, sem, 'k.', 'LineWidth',1.4);
    set(gca,'XTick',x,'XTickLabel',cellstr(groups));
    title(strrep(char(m),'_','\_'),'Interpreter','none');
    ylabel('Mean ± SEM');
    for gi=1:numel(groups)
        text(x(gi), vals(gi)+max(0.01,sem(gi)*1.2), sprintf('%.3f\nN=%d', vals(gi), round(ns(gi))), 'HorizontalAlignment','center', 'FontSize',8, 'Color',[0.2 0.2 0.2]);
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
bar(ax, x, vals, 0.62, 'FaceColor',[0.79 0.89 0.81], 'EdgeColor','none');
errorbar(ax, x, vals, sem, 'k.', 'LineWidth',1.4);
yline(ax, 0, '--', 'Color',[0.4 0.4 0.4]);
set(ax,'XTick',x,'XTickLabel',cellstr(M));
ylabel(ax,'Mean ± SEM');
title(ax, sprintf('Overall recovery summary [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
for i=1:numel(x)
    text(ax, x(i), vals(i)+sign(vals(i)+eps)*max(0.01,sem(i)*1.2), sprintf('%.3f\nN=%d', vals(i), round(T.N(i))), 'HorizontalAlignment','center', 'FontSize',8, 'Color',[0.2 0.2 0.2]);
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
    bar(x(i), vals(i), 0.62, 'FaceColor', colors(min(i,2),:), 'EdgeColor','none');
end
errorbar(ax, x, vals, sem, 'k.', 'LineWidth',1.4);
yline(ax, 0, '--', 'Color',[0.4 0.4 0.4]);
set(ax,'XTick',x,'XTickLabel',cellstr(groups));
ylabel(ax,'Mean delta_O_alpha ± SEM');
title(ax, sprintf('Experience recovery summary [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
for i=1:numel(x)
    text(ax, x(i), vals(i)+sign(vals(i)+eps)*max(0.01,sem(i)*1.2), sprintf('%.3f\nN=%d', vals(i), round(T.N(i))), 'HorizontalAlignment','center', 'FontSize',8, 'Color',[0.2 0.2 0.2]);
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
title(tl, sprintf('Experience descriptive summary by WWR × Complexity [%s]', tag), 'Interpreter','none', 'FontWeight','normal');
metrics = unique(string(T.metric), 'stable');
for mi=1:min(numel(metrics),4)
    m = metrics(mi);
    nexttile;
    hold on;
    style_axes(gca);
    X = T(string(T.metric)==m,:);
    plot_factor_lines(X, true);
    title(strrep(char(m),'_','\_'),'Interpreter','none');
end
pipeline.export_figure_png(fig, fullfile(fp_fig, sprintf('experience_factor_grid_%s.png', tag)), get_dpi(cfg));
try; close(fig); catch; end
end

function plot_factor_lines(X, hasGroup)
colors = [0.18 0.49 0.72; 0.88 0.47 0.18; 0.20 0.62 0.52; 0.62 0.42 0.78];
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
            [ww,o] = sort(str2double(string(Ti.WWR))); Ti = Ti(o,:);
            ls = '-'; if contains(lower(cx),'1') || contains(lower(cx),'high'), ls='--'; end
            errorbar(ww, Ti.mean, Ti.sem, 'o', 'LineWidth',1.5, 'Color', colors(min(ci,size(colors,1)),:), ...
                'MarkerFaceColor', colors(min(ci,size(colors,1)),:), 'DisplayName', sprintf('%s | %s', g, cx));
            plot(ww, Ti.mean, ls, 'LineWidth',1.7, 'Color', colors(min(ci,size(colors,1)),:), 'HandleVisibility','off');
            annotate_points(ww, Ti.mean, Ti.sem, Ti.N);
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
        [ww,o] = sort(str2double(string(Ti.WWR))); Ti = Ti(o,:);
        ls = '-'; if contains(lower(cx),'1') || contains(lower(cx),'high'), ls='--'; end
        errorbar(ww, Ti.mean, Ti.sem, 'o', 'LineWidth',1.5, 'Color', colors(cxi,:), ...
            'MarkerFaceColor', colors(cxi,:), 'DisplayName', sprintf('Complexity %s', cx));
        plot(ww, Ti.mean, ls, 'LineWidth',1.7, 'Color', colors(cxi,:), 'HandleVisibility','off');
        annotate_points(ww, Ti.mean, Ti.sem, Ti.N);
    end
    legend('Location','best');
end
xlabel('WWR');
ylabel('Mean ± SEM');
end

function annotate_points(x, y, se, n)
for i=1:numel(x)
    if ~isfinite(y(i)), continue; end
    yy = y(i) + max(0.01, abs(se(i))*1.2 + 0.01*max(1,abs(y(i))));
    text(x(i), yy, sprintf('%.3f\nN=%d', y(i), round(n(i))), 'HorizontalAlignment','center', 'FontSize',8, 'Color',[0.2 0.2 0.2]);
end
end

function style_axes(ax)
set(ax, 'Box','off', 'LineWidth',0.8, 'FontName','Arial', 'FontSize',11, 'Color','w');
grid(ax, 'on');
ax.GridAlpha = 0.12;
ax.GridColor = [0 0 0];
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
