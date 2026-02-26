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
fp_root = fullfile(fp_sum, 'analysis-2', 'task2_C1W45_block_diff');

branches = ["all","experience","sportfreq"];
fp_tbl = struct(); fp_rep = struct(); fp_fig = struct();
for bi=1:numel(branches)
    b = char(branches(bi));
    fp_tbl.(b) = fullfile(fp_root, 'tables', tag, b);
    fp_rep.(b) = fullfile(fp_root, 'reports', tag, b);
    fp_fig.(b) = fullfile(fp_root, 'figures', tag, b);
    if ~exist(fp_tbl.(b),'dir'); mkdir(fp_tbl.(b)); end
    if ~exist(fp_rep.(b),'dir'); mkdir(fp_rep.(b)); end
    if ~exist(fp_fig.(b),'dir'); mkdir(fp_fig.(b)); end
end

% required columns
req = {'subject_id','cond','block_id'};
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
try
    T.cond = string(T.cond);
end
T = T(lower(strtrim(string(T.cond)))=="view", :);

% select scene
sceneMask = true(height(T),1);
usedSceneKey = "scene_name";
if ismember('scene_name', T.Properties.VariableNames)
    sceneMask = (strtrim(string(T.scene_name))==string(scene_name));
else
    % fallback: if no scene_name, assume scene_id==1 corresponds to WWR45_C1 (as in current logs)
    if ismember('scene_id', T.Properties.VariableNames)
        sceneMask = (double(T.scene_id)==1);
        usedSceneKey = "scene_id==1";
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
writetable(Stats_all, fp_stats_all);
try
    plot_all_figures(S, metrics, fp_fig.all, tag, scene_name);
catch ME
    fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting all figures: %s\n', ME.message);
end
fp_md_all = fullfile(fp_rep.all, sprintf('scene_blockdiff_report_%s.md', tag));
write_report_md(fp_md_all, tag, metrics, scene_name, usedSceneKey, fp_sub_all, fp_stats_all, Stats_all, fp_fig.all, "all");

% --- branch: Experience groups ---
if ismember('Experience', S.Properties.VariableNames) && any(strlength(strtrim(string(S.Experience)))>0)
    fp_sub_ex = fullfile(fp_tbl.experience, sprintf('scene_blockdiff_subjectlevel_%s.csv', tag));
    writetable(S, fp_sub_ex);
    Stats_ex = pipeline.compute_blockdiff_stats(S, metrics, "group", "Experience");
    fp_stats_ex = fullfile(fp_tbl.experience, sprintf('scene_blockdiff_stats_%s.csv', tag));
    writetable(Stats_ex, fp_stats_ex);
    try
        plot_group_figures(S, metrics, fp_fig.experience, tag, scene_name, "Experience");
    catch ME
        fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting Experience figures: %s\n', ME.message);
    end
    fp_md_ex = fullfile(fp_rep.experience, sprintf('scene_blockdiff_report_%s.md', tag));
    write_report_md(fp_md_ex, tag, metrics, scene_name, usedSceneKey, fp_sub_ex, fp_stats_ex, Stats_ex, fp_fig.experience, "experience");
else
    Stats_ex = table(); fp_sub_ex = ""; fp_stats_ex = ""; fp_md_ex = "";
end

% --- branch: SportFreq groups ---
if ismember('SportFreq', S.Properties.VariableNames) && any(strlength(strtrim(string(S.SportFreq)))>0)
    fp_sub_sf = fullfile(fp_tbl.sportfreq, sprintf('scene_blockdiff_subjectlevel_%s.csv', tag));
    writetable(S, fp_sub_sf);
    Stats_sf = pipeline.compute_blockdiff_stats(S, metrics, "group", "SportFreq");
    fp_stats_sf = fullfile(fp_tbl.sportfreq, sprintf('scene_blockdiff_stats_%s.csv', tag));
    writetable(Stats_sf, fp_stats_sf);
    try
        plot_group_figures(S, metrics, fp_fig.sportfreq, tag, scene_name, "SportFreq");
    catch ME
        fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting SportFreq figures: %s\n', ME.message);
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
