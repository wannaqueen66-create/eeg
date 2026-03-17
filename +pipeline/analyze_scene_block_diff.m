function out = analyze_scene_block_diff(AllScene, fp_sum, cfg, tag, scene_name)
%ANALYZE_SCENE_BLOCK_DIFF Compare one specific scene between Block1 vs Block2.
%
% Default output root: <summary>/analysis-2/task2_C1W45_block_diff/
% (kept generic: scene_name is parameter; folder name fixed per user's request)
%
% For each subject and metric (default: O_theta, F_theta, O_alpha, O_beta),
% compute: diff = Block2 - Block1 for the target scene.
%
% Writes (per tag raw/qc) with 3 branches:
%   all/         (all subjects)
%   experience/  (Experience High/Low)
%   sportfreq/   (SportFreq High/Low)
% Each branch writes:
%   tables/<branch>/scene_blockdiff_subjectlevel_<tag>.csv
%   tables/<branch>/scene_blockdiff_stats_<tag>.csv
%   reports/<branch>/scene_blockdiff_report_<tag>.md
%   figures/<branch>/*.png (paper-friendly paired plots)

if nargin < 4 || isempty(tag)
    tag = 'raw';
end
if nargin < 3
    cfg = struct();
end
if nargin < 5 || strlength(string(scene_name))==0
    scene_name = "WWR45_C1";
end

out = struct();

% output dirs
[fp_root_base, fp_tbl_base, fp_fig_base, fp_rep_base] = pipeline.get_analysis_task_subdirs(fp_sum, 'task2_C1W45_block_diff', tag);

branches = ["all","experience","sportfreq"];
fp_tbl = struct(); fp_rep = struct(); fp_fig = struct();
for bi=1:numel(branches)
    b = char(branches(bi));
    fp_tbl.(b) = fullfile(fp_tbl_base, b);
    fp_rep.(b) = fullfile(fp_rep_base, b);
    fp_fig.(b) = fullfile(fp_fig_base, b);
    if ~exist(fp_tbl.(b),'dir'); mkdir(fp_tbl.(b)); end
    if ~exist(fp_rep.(b),'dir'); mkdir(fp_rep.(b)); end
    if ~exist(fp_fig.(b),'dir'); mkdir(fp_fig.(b)); end
end

% required columns
req = {'subject_id','block_id'};
for i=1:numel(req)
    if ~ismember(req{i}, AllScene.Properties.VariableNames)
        warning('analyze_scene_block_diff: missing required column %s. Skipping.', req{i});
        return;
    end
end

metrics = {"O_theta","F_theta","O_alpha","O_beta"};
try
    if isfield(cfg,'scene_blockdiff_metrics') && ~isempty(cfg.scene_blockdiff_metrics)
        metrics = string(cfg.scene_blockdiff_metrics);
    end
catch
end

T = AllScene;
% scene_level tables are view-only; `cond` may not exist here

% select scene
sceneMask = true(height(T),1);
usedSceneKey = "scene_name";
if ismember('scene_name', T.Properties.VariableNames)
    % case-insensitive robust match
    sceneMask = strcmpi(strtrim(string(T.scene_name)), strtrim(string(scene_name)));
    usedSceneKey = "scene_name (case-insensitive)";
else
    % Without scene_name, force explicit scene_id from config instead of hard-coded guess.
    if ismember('scene_id', T.Properties.VariableNames)
        scene_id_target = NaN;
        try
            if isfield(cfg,'task2_scene_id') && ~isempty(cfg.task2_scene_id)
                scene_id_target = double(cfg.task2_scene_id);
            end
        catch
        end
        if ~isfinite(scene_id_target)
            warning(['analyze_scene_block_diff: scene_name column missing, and cfg.task2_scene_id not provided. ', ...
                'Refusing implicit fallback to scene_id==1 to avoid mislabeling target scene.']);
            return;
        end
        sceneMask = (double(T.scene_id)==scene_id_target);
        usedSceneKey = sprintf('scene_id==%d (cfg.task2_scene_id)', round(scene_id_target));
    else
        warning('analyze_scene_block_diff: missing scene_name and scene_id; cannot isolate target scene.');
        return;
    end
end
T = T(sceneMask,:);

if isempty(T) || height(T)==0
    warning('analyze_scene_block_diff: no rows found for scene %s (%s).', string(scene_name), usedSceneKey);
    return;
end

% normalize factors if present
if ismember('Experience', T.Properties.VariableNames)
    T.Experience = normalize_high_low(T.Experience);
end
if ismember('SportFreq', T.Properties.VariableNames)
    T.SportFreq = normalize_high_low(T.SportFreq);
end

subs = unique(string(T.subject_id), 'stable');
rows = {};
for si=1:numel(subs)
    sid = subs(si);
    Ts = T(string(T.subject_id)==sid,:);

    iB1 = (double(Ts.block_id)==1);
    iB2 = (double(Ts.block_id)==2);
    if ~any(iB1) || ~any(iB2)
        continue;
    end

    ex = ""; sf = "";
    try
        if ismember('Experience', Ts.Properties.VariableNames)
            ex = string(Ts.Experience(find(iB1,1,'first')));
        end
    catch
    end
    try
        if ismember('SportFreq', Ts.Properties.VariableNames)
            sf = string(Ts.SportFreq(find(iB1,1,'first')));
        end
    catch
    end

    for mi=1:numel(metrics)
        m = string(metrics(mi));
        if ~ismember(m, Ts.Properties.VariableNames)
            continue;
        end
        y1 = double(Ts.(m)(find(iB1,1,'first')));
        y2 = double(Ts.(m)(find(iB2,1,'first')));
        if isnan(y1) || isnan(y2)
            continue;
        end
        d = y2 - y1;
        rows(end+1,:) = {sid, m, y1, y2, d, ex, sf}; %#ok<AGROW>
    end
end

if isempty(rows)
    warning('analyze_scene_block_diff: no usable subject pairs found for %s.', string(scene_name));
    return;
end

S = cell2table(rows, 'VariableNames', {'subject_id','metric','block1','block2','diff','Experience','SportFreq'});
S.scene_name = repmat(string(scene_name), height(S), 1);

% --- branch: all subjects ---
fp_sub_all = fullfile(fp_tbl.all, sprintf('scene_blockdiff_subjectlevel_%s.csv', tag));
writetable(S, fp_sub_all);
Stats_all = pipeline.compute_blockdiff_stats(S, metrics, "all", "");
fp_stats_all = fullfile(fp_tbl.all, sprintf('scene_blockdiff_stats_%s.csv', tag));
% Holm across metrics within this branch/test type
try
    Stats_all.p_ttest_holm = pipeline.holm_stepdown(double(Stats_all.p_ttest));
catch
    Stats_all.p_ttest_holm = nan(height(Stats_all),1);
end
try
    Stats_all.p_signrank_holm = pipeline.holm_stepdown(double(Stats_all.p_signrank));
catch
    Stats_all.p_signrank_holm = nan(height(Stats_all),1);
end
writetable(Stats_all, fp_stats_all);
try
    plot_all_figures(S, metrics, fp_fig.all, tag, scene_name);
catch ME
    fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting all figures: %s\n', ME.message);
end
try
    if ~isempty(Stats_all) && height(Stats_all)>0
        fp_ov_all = fullfile(fp_fig.all, pipeline.sanitize_filename(sprintf('scene_blockdiff_overview_%s_all_%s.png', tag, scene_name)));
        plot_branch_overview(Stats_all, metrics, tag, scene_name, "all", fp_ov_all);
    end
catch ME
    fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting all overview: %s\n', ME.message);
end
fp_md_all = fullfile(fp_rep.all, sprintf('scene_blockdiff_report_%s.md', tag));
write_report_md(fp_md_all, tag, metrics, scene_name, usedSceneKey, fp_sub_all, fp_stats_all, Stats_all, fp_fig.all, "all");

% --- branch: Experience groups ---
if ismember('Experience', S.Properties.VariableNames) && any(strlength(strtrim(string(S.Experience)))>0)
    fp_sub_ex = fullfile(fp_tbl.experience, sprintf('scene_blockdiff_subjectlevel_%s.csv', tag));
    writetable(S, fp_sub_ex);
    Stats_ex = pipeline.compute_blockdiff_stats(S, metrics, "group", "Experience");
    fp_stats_ex = fullfile(fp_tbl.experience, sprintf('scene_blockdiff_stats_%s.csv', tag));
    % Holm across metrics within each group (Experience) and test type
try
    Stats_ex.p_ttest_holm = nan(height(Stats_ex),1);
    Stats_ex.p_signrank_holm = nan(height(Stats_ex),1);
    for g = unique(string(Stats_ex.Group), 'stable')'
        idx = string(Stats_ex.Group)==g;
        Stats_ex.p_ttest_holm(idx) = pipeline.holm_stepdown(double(Stats_ex.p_ttest(idx)));
        Stats_ex.p_signrank_holm(idx) = pipeline.holm_stepdown(double(Stats_ex.p_signrank(idx)));
    end
catch
end
writetable(Stats_ex, fp_stats_ex);
    try
        plot_group_figures(S, metrics, fp_fig.experience, tag, scene_name, "Experience");
    catch ME
        fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting Experience figures: %s\n', ME.message);
    end
    try
        if ~isempty(Stats_ex) && height(Stats_ex)>0
            fp_ov_ex = fullfile(fp_fig.experience, pipeline.sanitize_filename(sprintf('scene_blockdiff_overview_%s_experience_%s.png', tag, scene_name)));
            plot_branch_overview(Stats_ex, metrics, tag, scene_name, "experience", fp_ov_ex);
        end
    catch ME
        fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting Experience overview: %s\n', ME.message);
    end
    fp_md_ex = fullfile(fp_rep.experience, sprintf('scene_blockdiff_report_%s.md', tag));
    write_report_md(fp_md_ex, tag, metrics, scene_name, usedSceneKey, fp_sub_ex, fp_stats_ex, Stats_ex, fp_fig.experience, "experience");
else
    Stats_ex = table(); fp_sub_ex = ""; fp_stats_ex = ""; fp_md_ex = "";
end

% --- branch: SportFreq groups ---
if ~isfield(fp_tbl, 'sportfreq')
    fp_tbl.sportfreq = fullfile(fp_tbl_base, 'sportfreq');
    if ~exist(fp_tbl.sportfreq,'dir'); mkdir(fp_tbl.sportfreq); end
end
if ~isfield(fp_rep, 'sportfreq')
    fp_rep.sportfreq = fullfile(fp_rep_base, 'sportfreq');
    if ~exist(fp_rep.sportfreq,'dir'); mkdir(fp_rep.sportfreq); end
end
if ~isfield(fp_fig, 'sportfreq')
    fp_fig.sportfreq = fullfile(fp_fig_base, 'sportfreq');
    if ~exist(fp_fig.sportfreq,'dir'); mkdir(fp_fig.sportfreq); end
end
if ismember('SportFreq', S.Properties.VariableNames) && any(strlength(strtrim(string(S.SportFreq)))>0)
    fp_sub_sf = fullfile(fp_tbl.sportfreq, sprintf('scene_blockdiff_subjectlevel_%s.csv', tag));
    writetable(S, fp_sub_sf);
    Stats_sf = pipeline.compute_blockdiff_stats(S, metrics, "group", "SportFreq");
    fp_stats_sf = fullfile(fp_tbl.sportfreq, sprintf('scene_blockdiff_stats_%s.csv', tag));
    % Holm across metrics within each group (SportFreq) and test type
try
    Stats_sf.p_ttest_holm = nan(height(Stats_sf),1);
    Stats_sf.p_signrank_holm = nan(height(Stats_sf),1);
    for g = unique(string(Stats_sf.Group), 'stable')'
        idx = string(Stats_sf.Group)==g;
        Stats_sf.p_ttest_holm(idx) = pipeline.holm_stepdown(double(Stats_sf.p_ttest(idx)));
        Stats_sf.p_signrank_holm(idx) = pipeline.holm_stepdown(double(Stats_sf.p_signrank(idx)));
    end
catch
end
writetable(Stats_sf, fp_stats_sf);
    try
        plot_group_figures(S, metrics, fp_fig.sportfreq, tag, scene_name, "SportFreq");
    catch ME
        fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting SportFreq figures: %s\n', ME.message);
    end
    try
        if ~isempty(Stats_sf) && height(Stats_sf)>0
            fp_ov_sf = fullfile(fp_fig.sportfreq, pipeline.sanitize_filename(sprintf('scene_blockdiff_overview_%s_sportfreq_%s.png', tag, scene_name)));
            plot_branch_overview(Stats_sf, metrics, tag, scene_name, "sportfreq", fp_ov_sf);
        end
    catch ME
        fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting SportFreq overview: %s\n', ME.message);
    end
    fp_md_sf = fullfile(fp_rep.sportfreq, sprintf('scene_blockdiff_report_%s.md', tag));
    write_report_md(fp_md_sf, tag, metrics, scene_name, usedSceneKey, fp_sub_sf, fp_stats_sf, Stats_sf, fp_fig.sportfreq, "sportfreq");
else
    Stats_sf = table(); fp_sub_sf = ""; fp_stats_sf = ""; fp_md_sf = "";
end

out.subjectlevel_csv = fp_sub_all;
out.stats_csv = fp_stats_all;
out.report_md = fp_md_all;
out.fig_dir = fp_fig.all;
out.branch_all = struct('subjectlevel_csv',fp_sub_all,'stats_csv',fp_stats_all,'report_md',fp_md_all,'fig_dir',fp_fig.all);
out.branch_experience = struct('subjectlevel_csv',fp_sub_ex,'stats_csv',fp_stats_ex,'report_md',fp_md_ex,'fig_dir',fp_fig.experience);
out.branch_sportfreq = struct('subjectlevel_csv',fp_sub_sf,'stats_csv',fp_stats_sf,'report_md',fp_md_sf,'fig_dir',fp_fig.sportfreq);
end

function y = normalize_high_low(x)
s = lower(strtrim(string(x)));
s(ismember(s,["high","h","1","true","yes","y","高"])) = "high";
s(ismember(s,["low","l","0","false","no","n","低"])) = "low";
y = repmat("", numel(s), 1);
y(s=="high") = "High";
y(s=="low")  = "Low";
end

function plot_all_figures(S, metrics, fp_fig, tag, scene_name)
if nargin<5; scene_name=""; end
for mi=1:numel(metrics)
    m = string(metrics(mi));
    X = S(S.metric==m,:);
    if height(X) < 3
        continue;
    end
    ttl = sprintf('%s | %s | All subjects | Block2 - Block1 [%s]', scene_name, m, tag);
    fp = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('scene_blockdiff_%s_all_%s_%s.png', tag, m, scene_name)));
    pipeline.plot_paired_summary(double(X.block1), double(X.block2), 'Block1','Block2', ttl, char(m), fp);
end
end

function plot_group_figures(S, metrics, fp_fig, tag, scene_name, gcol)
if nargin<5; scene_name=""; end
if nargin<6; gcol="Experience"; end

for g = ["Low","High"]
    for mi=1:numel(metrics)
        m = string(metrics(mi));
        X = S(S.metric==m & string(S.(gcol))==g, :);
        if height(X) < 3
            continue;
        end
        ttl = sprintf('%s | %s | %s=%s | Block2 - Block1 [%s]', scene_name, m, gcol, g, tag);
        fp = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('scene_blockdiff_%s_%s_%s_%s_%s.png', tag, lower(gcol), lower(g), m, scene_name)));
        pipeline.plot_paired_summary(double(X.block1), double(X.block2), 'Block1','Block2', ttl, char(m), fp);
    end
end
end

function plot_branch_overview(Stats, metrics, tag, scene_name, branchName, fp_png)
% Paper-friendly task2 overview per branch
% - rows: all (one row) or group rows (Low/High)
% - cols: metrics
% - cell color: mean(diff=Block2-Block1)
% - text: N, p_ttest, p_signrank, dz

set(0,'DefaultFigureVisible','off');

M = string(metrics(:));
nC = numel(M);

if isempty(Stats)
    return;
end

% Build row labels from available GroupType/Group in stats
if all(ismember({'GroupType','Group'}, Stats.Properties.VariableNames))
    rowKeys = unique(string(Stats.GroupType)+"|"+string(Stats.Group), 'stable');
else
    rowKeys = "all|All";
end
nR = numel(rowKeys);

Z = nan(nR,nC); P = nan(nR,nC); PSR = nan(nR,nC); DZ = nan(nR,nC); NN = nan(nR,nC);
for r=1:nR
    parts = split(rowKeys(r),'|');
    gt = parts(1); g = parts(2);
    for c=1:nC
        m = M(c);
        idx = string(Stats.metric)==m;
        if ismember('GroupType', Stats.Properties.VariableNames)
            idx = idx & string(Stats.GroupType)==gt;
        end
        if ismember('Group', Stats.Properties.VariableNames)
            idx = idx & string(Stats.Group)==g;
        end
        if any(idx)
            i = find(idx,1,'first');
            if ismember('mean_diff', Stats.Properties.VariableNames), Z(r,c)=double(Stats.mean_diff(i)); end
            if ismember('p_ttest', Stats.Properties.VariableNames), P(r,c)=double(Stats.p_ttest(i)); end
            if ismember('p_signrank', Stats.Properties.VariableNames), PSR(r,c)=double(Stats.p_signrank(i)); end
            if ismember('cohen_dz', Stats.Properties.VariableNames), DZ(r,c)=double(Stats.cohen_dz(i)); end
            if ismember('N', Stats.Properties.VariableNames), NN(r,c)=double(Stats.N(i)); end
        end
    end
end

fig = figure('Color','w','Position',[90 90 1200 480]);
ax = axes(fig); hold(ax,'on');
imagesc(ax, Z); axis(ax,'tight'); set(ax,'YDir','normal');
colormap(ax, interp1([0 0.5 1], [0.72 0.34 0.26; 0.96 0.96 0.96; 0.31 0.47 0.67], linspace(0,1,256)));
cb = colorbar(ax); cb.Label.String = 'mean diff = Block2 - Block1';
mx = max(abs(Z(:)),[],'omitnan');
if isempty(mx) || ~isfinite(mx) || mx==0, mx = 0.01; end
caxis(ax,[-mx mx]);

set(ax,'XTick',1:nC,'XTickLabel',cellstr(M));
set(ax,'YTick',1:nR,'YTickLabel',cellstr(replace(rowKeys,'|',' | ')));
xtickangle(ax,20);

title(ax, sprintf('Task2 scene block-diff overview | %s | %s [%s]', scene_name, branchName, tag), 'Interpreter','none');

for r=1:nR
    for c=1:nC
        if isnan(Z(r,c)); continue; end
        star = '';
        if ~isnan(P(r,c))
            if P(r,c) < 0.001, star='***';
            elseif P(r,c) < 0.01, star='**';
            elseif P(r,c) < 0.05, star='*';
            end
        end
        % display Holm-adjusted p (ttest) if available
        ph = P(r,c);
        try
            if ismember('p_ttest_holm', Stats.Properties.VariableNames)
                % locate the matching row again
                idx2 = string(Stats.metric)==m;
                if ismember('GroupType', Stats.Properties.VariableNames)
                    idx2 = idx2 & string(Stats.GroupType)==gt;
                end
                if ismember('Group', Stats.Properties.VariableNames)
                    idx2 = idx2 & string(Stats.Group)==g;
                end
                if any(idx2)
                    ph = double(Stats.p_ttest_holm(find(idx2,1,'first')));
                end
            end
        catch
        end
        txt = sprintf('N=%d\nΔ=%.3g%s\np_holm=%.3g\nsr=%.3g\ndz=%.2f', ...
            round(NN(r,c)), Z(r,c), star, ph, PSR(r,c), DZ(r,c));
        text(ax,c,r,txt,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
    end
end

pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function write_report_md(fp_md, tag, metrics, scene_name, usedSceneKey, fp_sub, fp_stats, Stats, fp_fig, branchName)
lines = {};
if nargin<10; branchName=""; end
lines{end+1} = sprintf('# Scene block difference (%s) [%s]', tag, string(branchName));
lines{end+1} = '';
lines{end+1} = sprintf('Target scene: **%s** (selected by %s)', string(scene_name), string(usedSceneKey));
lines{end+1} = '';
lines{end+1} = 'Definition:';
lines{end+1} = '- Use **view** segments for the target scene only.';
lines{end+1} = '- For each subject and metric: **diff = Block2 − Block1**.';
lines{end+1} = '';
lines{end+1} = 'Outputs:';
lines{end+1} = sprintf('- Subject-level: `%s`', fp_sub);
lines{end+1} = sprintf('- Stats: `%s`', fp_stats);
lines{end+1} = sprintf('- Figures: `%s`', fp_fig);
lines{end+1} = sprintf('- Overview figure: `%s`', fullfile(fp_fig, pipeline.sanitize_filename(sprintf('scene_blockdiff_overview_%s_%s_%s.png', tag, branchName, scene_name))));
lines{end+1} = '';
lines{end+1} = 'Metrics:';
lines{end+1} = sprintf('- %s', strjoin(string(metrics), ', '));
lines{end+1} = '';

if isempty(Stats) || height(Stats)==0
    lines{end+1} = '> No group stats computed (missing grouping columns or insufficient N).';
else
    lines{end+1} = 'Quick summary (paired t-test p-values):';
    lines{end+1} = '';
    for i=1:height(Stats)
        lines{end+1} = sprintf('- %s %s | %s: N=%d, mean(diff)=%.4g, t(%d)=%.3f, p=%.3g (signrank p=%.3g)', ...
            string(Stats.GroupType(i)), string(Stats.Group(i)), string(Stats.metric(i)), ...
            Stats.N(i), Stats.mean_diff(i), Stats.df(i), Stats.t(i), Stats.p_ttest(i), Stats.p_signrank(i));
    end
end

fid = fopen(fp_md,'w');
if fid<0
    warning('analyze_scene_block_diff: cannot write report %s', fp_md);
    return;
end
fwrite(fid, strjoin(string(lines), newline));
fclose(fid);
end
