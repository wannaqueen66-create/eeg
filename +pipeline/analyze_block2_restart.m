function out = analyze_block2_restart(AllScene, fp_sum, cfg, tag)
%ANALYZE_BLOCK2_RESTART Block2 restart effect: compare Block2 first scene vs mean of remaining 5.
%
% Uses pipeline.plot_paired_summary for figures.
%
% Output root (requested by user): <summary>/analysis-2/task1_block2_restart/
% Writes (per tag: raw/qc):
%   tables/<tag>/block2_restart_subjectlevel_<tag>.csv
%   tables/<tag>/block2_restart_stats_<tag>.csv
%   reports/<tag>/block2_restart_report_<tag>.md
%   figures/<tag>/*.png (paper-friendly paired plots)
%
% Requires columns in AllScene:
%   subject_id, block_id, cycle_in_block
% Note: summarize_bandpower_outputs merges per-subject *_scene_level.csv which are view-only and do NOT include `cond`.
%   metrics: O_theta, F_theta, O_alpha, O_beta (configurable)
% Optional grouping columns:
%   Experience, SportFreq  (values: High/Low or 高/低 or 1/0)

if nargin < 4 || isempty(tag)
    tag = 'raw';
end
if nargin < 3
    cfg = struct();
end

out = struct();

% --- output dirs ---
[fp_root, fp_tbl, fp_fig, fp_rep] = pipeline.get_analysis_task_subdirs(fp_sum, 'task1_block2_restart', tag);

% --- required columns ---
req = {'subject_id','block_id','cycle_in_block'};
for i=1:numel(req)
    if ~ismember(req{i}, AllScene.Properties.VariableNames)
        warning('analyze_block2_restart: missing required column %s. Skipping.', req{i});
        return;
    end
end

% metrics list (default requested)
metrics = {"O_theta","F_theta","O_alpha","O_beta"};
try
    if isfield(cfg,'block2_restart_metrics') && ~isempty(cfg.block2_restart_metrics)
        metrics = string(cfg.block2_restart_metrics);
    end
catch
end

% filter Block2 (scene_level tables are already view-only)
T = AllScene;
T = T(double(T.block_id)==2, :);

if isempty(T) || height(T)==0
    warning('analyze_block2_restart: no Block2 view rows found.');
    return;
end

% normalize grouping columns (if present)
if ismember('Experience', T.Properties.VariableNames)
    T.Experience = normalize_high_low(T.Experience);
elseif ismember('ExperienceGroup', T.Properties.VariableNames)
    T.Experience = normalize_high_low(T.ExperienceGroup);
end
if ismember('SportFreq', T.Properties.VariableNames)
    T.SportFreq = normalize_high_low(T.SportFreq);
elseif ismember('SportFreqGroup', T.Properties.VariableNames)
    T.SportFreq = normalize_high_low(T.SportFreqGroup);
end

% subject-level wide to per-metric long
subs = unique(string(T.subject_id), 'stable');

rows = {};
for si=1:numel(subs)
    sid = subs(si);
    Ts = T(string(T.subject_id)==sid, :);

    % First scene in Block2
    iFirst = (double(Ts.cycle_in_block)==1);
    iRest  = (double(Ts.cycle_in_block)>=2 & double(Ts.cycle_in_block)<=6);

    if ~any(iFirst) || ~any(iRest)
        continue;
    end

    % group labels (may be missing)
    ex = ""; sf = "";
    try
        if ismember('Experience', Ts.Properties.VariableNames)
            ex = string(Ts.Experience(find(iFirst,1,'first')));
        end
    catch
    end
    try
        if ismember('SportFreq', Ts.Properties.VariableNames)
            sf = string(Ts.SportFreq(find(iFirst,1,'first')));
        end
    catch
    end

    for mi=1:numel(metrics)
        m = string(metrics(mi));
        if ~ismember(m, Ts.Properties.VariableNames)
            continue;
        end
        yFirst = Ts.(m)(iFirst);
        yRest  = Ts.(m)(iRest);
        yFirst = double(yFirst(1));
        yRestM = mean(double(yRest), 'omitnan');
        if isnan(yFirst) || isnan(yRestM)
            continue;
        end
        diffv = yFirst - yRestM;

        rows(end+1,:) = {sid, m, yFirst, yRestM, diffv, ex, sf}; %#ok<AGROW>
    end
end

if isempty(rows)
    warning('analyze_block2_restart: no usable subject-level pairs found.');
    return;
end

S = cell2table(rows, 'VariableNames', {'subject_id','metric','block2_first','block2_rest_mean','diff','Experience','SportFreq'});

% write subject-level table
fp_sub = fullfile(fp_tbl, sprintf('block2_restart_subjectlevel_%s.csv', tag));
writetable(S, fp_sub);

% stats per groupType x group x metric
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

            % paired t-test first vs rest
            try
                [~, p, ~, st] = ttest(double(X.block2_first), double(X.block2_rest_mean));
                tstat = st.tstat;
                df = st.df;
            catch
                p = NaN; tstat = NaN; df = NaN;
            end

            % signrank as robust alternative
            try
                p_sr = signrank(double(X.block2_first), double(X.block2_rest_mean));
            catch
                p_sr = NaN;
            end

            statsRows(end+1,:) = {gcol, g, m, n, md, sem, tstat, df, p, dz, p_sr}; %#ok<AGROW>
        end
    end
end

if isempty(statsRows)
    Stats = table('Size',[0 11], 'VariableTypes', ...
        {'string','string','string','double','double','double','double','double','double','double','double'}, ...
        'VariableNames', {'GroupType','Group','metric','N','mean_diff','sem_diff','t','df','p_ttest','cohen_dz','p_signrank'});
else
    Stats = cell2table(statsRows, 'VariableNames', {
        'GroupType','Group','metric','N','mean_diff','sem_diff','t','df','p_ttest','cohen_dz','p_signrank'});
end

% Holm across metrics within each group (GroupType×Group) and test type
try
    Stats.p_ttest_holm = nan(height(Stats),1);
    Stats.p_signrank_holm = nan(height(Stats),1);
    keys = unique(string(Stats.GroupType)+"|"+string(Stats.Group), 'stable');
    for ki=1:numel(keys)
        parts = split(keys(ki),'|');
        gt = parts(1); g = parts(2);
        idx = string(Stats.GroupType)==gt & string(Stats.Group)==g;
        Stats.p_ttest_holm(idx) = pipeline.holm_stepdown(double(Stats.p_ttest(idx)));
        Stats.p_signrank_holm(idx) = pipeline.holm_stepdown(double(Stats.p_signrank(idx)));
    end
catch
end

fp_stats = fullfile(fp_tbl, sprintf('block2_restart_stats_%s.csv', tag));
writetable(Stats, fp_stats);

% paper-friendly figures
try
    plot_paired_figures(S, metrics, fp_fig, tag);
catch ME
    fprintf(2, '[WARN] analyze_block2_restart: failed plotting figures: %s\n', ME.message);
end

% overview figure (paper-friendly summary of mean_diff / p / dz)
try
    if ~isempty(Stats) && height(Stats)>0
        fp_over = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('block2_restart_overview_%s.png', tag)));
        plot_overview_figure(Stats, metrics, tag, fp_over);
    end
catch ME
    fprintf(2, '[WARN] analyze_block2_restart: failed overview figure: %s\n', ME.message);
end

% report markdown
fp_md = fullfile(fp_rep, sprintf('block2_restart_report_%s.md', tag));
write_report_md(fp_md, tag, metrics, fp_sub, fp_stats, Stats, fp_fig);

out.subjectlevel_csv = fp_sub;
out.stats_csv = fp_stats;
out.report_md = fp_md;
out.fig_dir = fp_fig;
end

function y = normalize_high_low(x)
% Map various encodings -> "High"/"Low" (otherwise empty)
s = lower(strtrim(string(x)));

% common variants
s(ismember(s,["high","h","1","true","yes","y","高"])) = "high";
s(ismember(s,["low","l","0","false","no","n","低"])) = "low";

% final
y = repmat("", numel(s), 1);
y(s=="high") = "High";
y(s=="low")  = "Low";
end

function plot_paired_figures(S, metrics, fp_fig, tag)
% Paper-friendly paired plots per metric, within each group type and group.

set(0,'DefaultFigureVisible','off');

if nargin<4; tag='raw'; end

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
            y1 = double(X.block2_first);
            y2 = double(X.block2_rest_mean);
            ttl = sprintf('Block2 restart | %s | %s=%s | %s [%s]', m, gcol, g, 'first - mean(rest)', tag);
            fp = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('block2_restart_%s_%s_%s_%s.png', tag, lower(gcol), lower(g), m)));
            pipeline.plot_paired_summary(y1, y2, 'B2-1', 'B2-2..6 mean', ttl, char(m), fp);
        end
    end
end
end

function plot_overview_figure(Stats, metrics, tag, fp_png)
% Paper-friendly task1 overview:
% - rows: GroupType x Group
% - cols: metrics
% - cell color: mean_diff
% - text: N, p_ttest, p_signrank, dz

set(0,'DefaultFigureVisible','off');

GT = ["Experience","SportFreq"];
G  = ["Low","High"];
rows = strings(0,1);
for gi=1:numel(GT)
    for gj=1:numel(G)
        rows(end+1,1) = GT(gi)+"|"+G(gj); %#ok<AGROW>
    end
end
nR = numel(rows);
M = string(metrics(:));
nC = numel(M);

Z = nan(nR,nC);
P = nan(nR,nC);
PSR = nan(nR,nC);
DZ = nan(nR,nC);
NN = nan(nR,nC);

for r=1:nR
    parts = split(rows(r),'|');
    gt = parts(1); g = parts(2);
    for c=1:nC
        m = M(c);
        idx = string(Stats.GroupType)==gt & string(Stats.Group)==g & string(Stats.metric)==m;
        if any(idx)
            i = find(idx,1,'first');
            Z(r,c) = double(Stats.mean_diff(i));
            P(r,c) = double(Stats.p_ttest(i));
            PSR(r,c)= double(Stats.p_signrank(i));
            DZ(r,c)= double(Stats.cohen_dz(i));
            NN(r,c)= double(Stats.N(i));
        end
    end
end

fig = figure('Color','w','Position',[80 80 1300 560]);
ax = axes(fig); hold(ax,'on');

imagesc(ax, Z);
axis(ax,'tight');
set(ax,'YDir','normal');
colormap(ax, interp1([0 0.5 1], [0.72 0.34 0.26; 0.96 0.96 0.96; 0.31 0.47 0.67], linspace(0,1,256)));
cb = colorbar(ax);
cb.Label.String = 'mean diff = B2-1 - mean(B2-2..6)';

% Symmetric color scaling around 0 for sign readability
mx = max(abs(Z(:)),[],'omitnan');
if isempty(mx) || ~isfinite(mx) || mx==0
    mx = 0.01;
end
caxis(ax, [-mx mx]);

set(ax,'XTick',1:nC,'XTickLabel',cellstr(M));
set(ax,'YTick',1:nR,'YTickLabel',cellstr(replace(rows, '|', ' | ')));
xtickangle(ax, 20);

title(ax, sprintf('Task1 Block2 restart overview [%s]', tag), 'Interpreter','none');

% annotate each cell
for r=1:nR
    for c=1:nC
        if isnan(Z(r,c)); continue; end
        star = '';
        if ~isnan(P(r,c))
            if P(r,c) < 0.001, star = '***';
            elseif P(r,c) < 0.01, star = '**';
            elseif P(r,c) < 0.05, star = '*';
            end
        end
        % display Holm-adjusted p (ttest) if available
        ph = P(r,c);
        try
            if ismember('p_ttest_holm', Stats.Properties.VariableNames)
                idx2 = string(Stats.GroupType)==gt & string(Stats.Group)==g & string(Stats.metric)==m;
                if any(idx2)
                    ph = double(Stats.p_ttest_holm(find(idx2,1,'first')));
                end
            end
        catch
        end
        txt = sprintf('N=%d\nΔ=%.3g%s\np_holm=%.3g\nsr=%.3g\ndz=%.2f', ...
            round(NN(r,c)), Z(r,c), star, ph, PSR(r,c), DZ(r,c));
        text(ax, c, r, txt, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'Color','k', 'FontSize',8, 'FontWeight','normal');
    end
end

grid(ax,'off');
pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function write_report_md(fp_md, tag, metrics, fp_sub, fp_stats, Stats, fp_fig)
lines = {};
lines{end+1} = sprintf('# Block2 restart analysis (%s)', tag);
lines{end+1} = '';
lines{end+1} = 'Definition:';
lines{end+1} = '- Use **Block2 view** segments only.';
lines{end+1} = '- For each subject and metric: **diff = Block2(first scene) − mean(Block2 scenes 2–6)**.';
lines{end+1} = '';
lines{end+1} = 'Outputs:';
lines{end+1} = sprintf('- Subject-level table: `%s`', fp_sub);
lines{end+1} = sprintf('- Stats table: `%s`', fp_stats);
try
    if nargin>=7 && strlength(string(fp_fig))>0
        lines{end+1} = sprintf('- Figures dir: `%s`', fp_fig);
        lines{end+1} = sprintf('- Overview figure: `%s`', fullfile(fp_fig, pipeline.sanitize_filename(sprintf('block2_restart_overview_%s.png', tag))));
    end
catch
end
lines{end+1} = '';
lines{end+1} = 'Metrics:';
lines{end+1} = sprintf('- %s', strjoin(string(metrics), ', '));
lines{end+1} = '';

if isempty(Stats) || height(Stats)==0
    lines{end+1} = '> No group stats computed (missing grouping columns or insufficient N).';
else
    lines{end+1} = 'Quick summary (paired t-test p-values):';
    lines{end+1} = '';
    % write a compact bullet list
    for i=1:height(Stats)
        lines{end+1} = sprintf('- %s %s | %s: N=%d, mean(diff)=%.4g, t(%d)=%.3f, p=%.3g (signrank p=%.3g)', ...
            string(Stats.GroupType(i)), string(Stats.Group(i)), string(Stats.metric(i)), ...
            Stats.N(i), Stats.mean_diff(i), Stats.df(i), Stats.t(i), Stats.p_ttest(i), Stats.p_signrank(i));
    end
end

txt = strjoin(string(lines), newline);
fid = fopen(fp_md,'w');
if fid<0
    warning('analyze_block2_restart: cannot write report %s', fp_md);
    return;
end
fwrite(fid, txt);
fclose(fid);
end
