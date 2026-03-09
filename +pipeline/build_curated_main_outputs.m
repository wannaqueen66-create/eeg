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
        end
    end
    if ~isempty(AllScene_qc)
        writetable(AllScene_qc, fullfile(fp_do_tbl, 'scene_level_overall_qc.csv'));
        Tsum = summarize_overall_scene(AllScene_qc, cfg);
        if ~isempty(Tsum)
            writetable(Tsum, fullfile(fp_do_tbl, 'overall_scene_metric_means_qc.csv'));
        end
    end
    if ~isempty(AllPairs)
        writetable(AllPairs, fullfile(fp_do_tbl, 'pairs_overall_raw.csv'));
    end
    if ~isempty(AllPairs_qc)
        writetable(AllPairs_qc, fullfile(fp_do_tbl, 'pairs_overall_qc.csv'));
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
        fprintf(fid, '- `pairs_overall_raw.csv`: full-sample recovery/pair descriptive table\n');
        fprintf(fid, '- `pairs_overall_qc.csv`: QC-filtered recovery/pair descriptive table\n');
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
        fprintf(fid, '- `pairs_experience_*.csv`: recovery tables with valid experience labels\n');
        fprintf(fid, '- `experience_recovery_means_*.csv`: recovery summaries by Experience group\n');
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

function build_inferential_experience_summary(fp_sum, fp_tbl, fp_fig, fp_rep)
rows = {};
% Task4 model1/model2/model3 and Task3 key ANOVA summaries for experience branch
cands = {
    fullfile(fp_sum,'analysis','task4_core_lmm_suite','raw','tables','factor_WWR','experience'), ...
    fullfile(fp_sum,'analysis','task4_core_lmm_suite','qc','tables','factor_WWR','experience'), ...
    fullfile(fp_sum,'analysis','task3_trialindex_lmm','raw','tables','experience'), ...
    fullfile(fp_sum,'analysis','task3_trialindex_lmm','qc','tables','experience') ...
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
end
fid = fopen(fullfile(fp_rep,'README.md'),'w');
if fid>0
    fprintf(fid, '# Inferential / Experience\n\n');
    fprintf(fid, 'This folder is the curated entry point for experience-group inferential outputs.\n');
    fprintf(fid, 'Use `experience_inferential_file_index.csv` as the fastest file-level map into the detailed analysis outputs.\n\n');
    fprintf(fid, 'Primary source tasks currently include:\n');
    fprintf(fid, '- task3_trialindex_lmm\n');
    fprintf(fid, '- task4_core_lmm_suite\n');
    fprintf(fid, '- task5_peakindex_invertedu\n');
    fprintf(fid, '- task6_coremetric_special\n');
    fprintf(fid, '- task7_individual_checks\n');
    fclose(fid);
end
end
