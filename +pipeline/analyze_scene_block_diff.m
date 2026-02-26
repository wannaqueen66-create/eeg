function out = analyze_scene_block_diff(AllScene, fp_sum, cfg, tag, scene_name)
%ANALYZE_SCENE_BLOCK_DIFF Compare one specific scene between Block1 vs Block2.
%
% Default output root: <summary>/analysis-2/task2_C1W45_block_diff/
% (kept generic: scene_name is parameter; folder name fixed per user's request)
%
% For each subject and metric (default: O_theta, F_theta, O_alpha, O_beta),
% compute: diff = Block2 - Block1 for the target scene.
%
% Writes (per tag raw/qc):
%   tables/<tag>/scene_blockdiff_subjectlevel_<tag>.csv
%   tables/<tag>/scene_blockdiff_stats_<tag>.csv
%   reports/<tag>/scene_blockdiff_report_<tag>.md
%   figures/<tag>/*.png (paper-friendly paired plots)

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
fp_tbl = fullfile(fp_root, 'tables', tag);
fp_rep = fullfile(fp_root, 'reports', tag);
fp_fig = fullfile(fp_root, 'figures', tag);
if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
if ~exist(fp_fig,'dir'); mkdir(fp_fig); end

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

fp_sub = fullfile(fp_tbl, sprintf('scene_blockdiff_subjectlevel_%s.csv', tag));
writetable(S, fp_sub);

% stats per group
statsRows = {};
groupTypes = ["Experience","SportFreq"];
for gi=1:numel(groupTypes)
    gcol = groupTypes(gi);
    if ~ismember(gcol, S.Properties.VariableNames)
        continue;
    end
    for g = ["Low","High"]
        for mi=1:numel(metrics)
            m = string(metrics(mi));
            X = S(S.metric==m & string(S.(gcol))==g, :);
            if height(X) < 3
                continue;
            end
            d = double(X.diff);
            n = sum(~isnan(d));
            if n < 3
                continue;
            end
            md = mean(d,'omitnan');
            sd = std(d,'omitnan');
            sem = sd / sqrt(n);
            dz = md / sd;

            try
                [~, p, ~, st] = ttest(double(X.block2), double(X.block1));
                tstat = st.tstat;
                df = st.df;
            catch
                p = NaN; tstat = NaN; df = NaN;
            end
            try
                p_sr = signrank(double(X.block2), double(X.block1));
            catch
                p_sr = NaN;
            end

            statsRows(end+1,:) = {gcol, g, m, n, md, sem, tstat, df, p, dz, p_sr}; %#ok<AGROW>
        end
    end
end

Stats = cell2table(statsRows, 'VariableNames', {
    'GroupType','Group','metric','N','mean_diff','sem_diff','t','df','p_ttest','cohen_dz','p_signrank'});

fp_stats = fullfile(fp_tbl, sprintf('scene_blockdiff_stats_%s.csv', tag));
writetable(Stats, fp_stats);

% figures
try
    plot_paired_figures(S, metrics, fp_fig, tag, scene_name);
catch ME
    fprintf(2,'[WARN] analyze_scene_block_diff: failed plotting figures: %s\n', ME.message);
end

% report
fp_md = fullfile(fp_rep, sprintf('scene_blockdiff_report_%s.md', tag));
write_report_md(fp_md, tag, metrics, scene_name, usedSceneKey, fp_sub, fp_stats, Stats, fp_fig);

out.subjectlevel_csv = fp_sub;
out.stats_csv = fp_stats;
out.report_md = fp_md;
out.fig_dir = fp_fig;
end

function y = normalize_high_low(x)
s = lower(strtrim(string(x)));
s(ismember(s,["high","h","1","true","yes","y","高"])) = "high";
s(ismember(s,["low","l","0","false","no","n","低"])) = "low";
y = repmat("", numel(s), 1);
y(s=="high") = "High";
y(s=="low")  = "Low";
end

function plot_paired_figures(S, metrics, fp_fig, tag, scene_name)
if nargin<5; scene_name = ""; end

set(0,'DefaultFigureVisible','off');

% pick group types if present
groupTypes = ["Experience","SportFreq"];
for gi=1:numel(groupTypes)
    gcol = groupTypes(gi);
    if ~ismember(gcol, S.Properties.VariableNames)
        continue;
    end
    for g = ["Low","High"]
        for mi=1:numel(metrics)
            m = string(metrics(mi));
            X = S(S.metric==m & string(S.(gcol))==g, :);
            if height(X) < 3
                continue;
            end

            y1 = double(X.block1);
            y2 = double(X.block2);
            n = numel(y1);

            fig = figure('Position',[100 100 780 420], 'Color','w');
            hold on;
            for i=1:n
                plot([1 2], [y1(i) y2(i)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            end
            scatter(ones(n,1)*1, y1, 30, 'filled', 'MarkerFaceColor',[0.2 0.5 0.9]);
            scatter(ones(n,1)*2, y2, 30, 'filled', 'MarkerFaceColor',[0.9 0.4 0.2]);

            mu1 = mean(y1,'omitnan'); mu2 = mean(y2,'omitnan');
            se1 = std(y1,'omitnan')/sqrt(n); se2 = std(y2,'omitnan')/sqrt(n);
            errorbar(1, mu1, se1, 'k', 'LineWidth', 2);
            errorbar(2, mu2, se2, 'k', 'LineWidth', 2);

            xlim([0.6 2.4]);
            set(gca,'XTick',[1 2],'XTickLabel',{'Block1','Block2'});
            grid on;
            title(sprintf('%s | %s | %s=%s | %s [%s]', scene_name, m, gcol, g, 'Block2 - Block1', tag), 'Interpreter','none');
            ylabel(char(m));

            fn = pipeline.sanitize_filename(sprintf('scene_blockdiff_%s_%s_%s_%s_%s.png', tag, lower(gcol), lower(g), m, scene_name));
            fp = fullfile(fp_fig, fn);
            try
                exportgraphics(fig, fp, 'Resolution', 300);
            catch
                saveas(fig, fp);
            end
            try; close(fig); catch; end
        end
    end
end
end

function write_report_md(fp_md, tag, metrics, scene_name, usedSceneKey, fp_sub, fp_stats, Stats, fp_fig)
lines = {};
lines{end+1} = sprintf('# Scene block difference (%s)', tag);
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
