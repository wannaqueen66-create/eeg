function out = analyze_peakindex_invertedu(AllScene, fp_sum, cfg, tag)
%ANALYZE_PEAKINDEX_INVERTEDU Task5: Inverted-U PeakIndex test.
%
% For each metric (default: O_theta, F_theta), and each grouping analysis
% (ExperienceGroup / SportFreqGroup), compute subject-level:
%   PeakIndex = EEG(45) - mean(EEG(15), EEG(75))
% within each Complexity level.
%
% Then fit:
%   PeakIndex ~ Complexity*Group + (1|Subject)
%
% Outputs under:
%   <summary>/analysis-2/task5_peakindex_invertedU/
%     tables/<tag>/{experience|sportfreq}/...
%     reports/<tag>/{experience|sportfreq}/...
%     figures/<tag>/{experience|sportfreq}/...

if nargin < 4 || isempty(tag)
    tag = 'raw';
end
if nargin < 3
    cfg = struct();
end

out = struct();

if exist('fitlme','file') ~= 2
    warning('analyze_peakindex_invertedu: fitlme not found (Stats toolbox missing). Skipping.');
    return;
end

req = {'subject_id','WWR','Complexity'};
for i=1:numel(req)
    if ~ismember(req{i}, AllScene.Properties.VariableNames)
        warning('analyze_peakindex_invertedu: missing required column %s. Skipping.', req{i});
        return;
    end
end

metrics = ["O_theta","F_theta"];
try
    if isfield(cfg,'task5_metrics') && ~isempty(cfg.task5_metrics)
        metrics = string(cfg.task5_metrics);
    end
catch
end

analyses = { ...
    struct('name','experience','gcol','ExperienceGroup'), ...
    struct('name','sportfreq','gcol','SportFreqGroup') ...
};

fp_root = fullfile(fp_sum, 'analysis-2', 'task5_peakindex_invertedU');

T0 = AllScene;
T0.subject_id = string(T0.subject_id);
T0.WWR = normalize_wwr(T0.WWR);
T0.Complexity = normalize_complexity(T0.Complexity);

for ai = 1:numel(analyses)
    A = analyses{ai};
    gcol = string(A.gcol);
    if ~ismember(gcol, string(T0.Properties.VariableNames))
        continue;
    end
    gcolc = char(gcol);

    SumRows = table();

    for mi = 1:numel(metrics)
        dv = metrics(mi);
        if ~ismember(dv, string(T0.Properties.VariableNames))
            continue;
        end
        dvc = char(dv);

        T = T0(:, {'subject_id','WWR','Complexity',gcolc,dvc});
        T.Properties.VariableNames = {'subject_id','WWR','Complexity','GroupRaw','EEGraw'};
        T.EEG = double(T.EEGraw);
        T.Group = normalize_high_low(T.GroupRaw);
        T = T(~isnan(T.EEG), :);
        T = T(strlength(string(T.WWR))>0 & strlength(string(T.Complexity))>0 & strlength(string(T.Group))>0, :);
        T = T(ismember(string(T.WWR), ["15","45","75"]), :);
        if height(T) < 30
            continue;
        end

        T.Subject = categorical(string(T.subject_id));
        T.WWR = categorical(string(T.WWR), {'15','45','75'});
        T.Complexity = categorical(string(T.Complexity), {'ComplexityLow','ComplexityHigh'});
        T.Group = categorical(string(T.Group), {'Low','High'});

        % subject-level mean EEG per WWR within (Subject, Complexity, Group)
        [G, sid, grp, cx, wwr] = findgroups(T.Subject, T.Group, T.Complexity, T.WWR);
        mu = splitapply(@(x) mean(x,'omitnan'), T.EEG, G);
        Awwr = table(sid, grp, cx, string(wwr), mu, 'VariableNames', {'Subject','Group','Complexity','WWR','EEG'});

        % Wide by WWR level
        K = unique(Awwr(:, {'Subject','Group','Complexity'}), 'rows');
        for lv = ["15","45","75"]
            Z = Awwr(Awwr.WWR==lv, {'Subject','Group','Complexity','EEG'});
            Z.Properties.VariableNames{'EEG'} = ['EEG_' char(lv)];
            K = outerjoin(K, Z, 'Keys', {'Subject','Group','Complexity'}, 'MergeKeys', true, 'Type','left');
        end

        if ~ismember('EEG_45', K.Properties.VariableNames)
            continue;
        end
        if ~ismember('EEG_15', K.Properties.VariableNames)
            K.EEG_15 = nan(height(K),1);
        end
        if ~ismember('EEG_75', K.Properties.VariableNames)
            K.EEG_75 = nan(height(K),1);
        end

        K.PeakIndex = K.EEG_45 - mean([K.EEG_15 K.EEG_75], 2, 'omitnan');
        K = K(~isnan(K.PeakIndex), :);
        if height(K) < 12
            continue;
        end

        % Output dirs
        fp_tbl = fullfile(fp_root, 'tables', tag, A.name);
        fp_rep = fullfile(fp_root, 'reports', tag, A.name);
        fp_fig = fullfile(fp_root, 'figures', tag, A.name);
        if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
        if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
        if ~exist(fp_fig,'dir'); mkdir(fp_fig); end

        % Save subject-level PeakIndex
        fp_subject = fullfile(fp_tbl, sprintf('peakindex_subjectlevel_%s_%s.csv', dv, tag));
        writetable(K, fp_subject);

        % LMM
        try
            lme = fitlme(K, 'PeakIndex ~ Complexity*Group + (1|Subject)');
        catch ME
            warning('analyze_peakindex_invertedu: fitlme failed (%s/%s/%s): %s', tag, A.name, dv, ME.message);
            continue;
        end

        % export tables
        try
            writetable(lme.Coefficients, fullfile(fp_tbl, sprintf('peakindex_lmm_fixed_effects_%s_%s.csv', dv, tag)));
        catch
        end
        try
            Aov = anova(lme,'DFMethod','Satterthwaite');
            writetable(Aov, fullfile(fp_tbl, sprintf('peakindex_lmm_anova_%s_%s.csv', dv, tag)));
        catch
        end

        % report
        try
            fp_md = fullfile(fp_rep, sprintf('peakindex_report_%s_%s.md', dv, tag));
            write_report(fp_md, dv, tag, A.name, K, lme);
        catch
        end

        % figure
        try
            fp_png = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('peakindex_%s_%s_%s.png', tag, A.name, dv)));
            plot_peakindex_figure(K, dv, tag, A.name, fp_png);
        catch
        end

        % collect overview stats row
        try
            [nsub, mu_all, se_all, p0, bG, pG, bC, pC, bCG, pCG] = extract_peakindex_summary(K, lme);
            rr = table(string(A.name), dv, nsub, mu_all, se_all, p0, bG, pG, bC, pC, bCG, pCG, ...
                'VariableNames', {'analysis','metric','n_subjects','mean_peakindex','sem_peakindex','p_vs_zero', ...
                                  'beta_group','p_group','beta_complexity','p_complexity','beta_cxg','p_cxg'});
            SumRows = [SumRows; rr]; %#ok<AGROW>
        catch
        end
    end

    % overview outputs for this grouping analysis
    try
        if ~isempty(SumRows) && height(SumRows)>0
            fp_tbl_over = fullfile(fp_root, 'tables', tag, A.name);
            fp_fig_over = fullfile(fp_root, 'figures', tag, A.name);
            fp_csv_over = fullfile(fp_tbl_over, sprintf('peakindex_summary_%s.csv', tag));
            writetable(SumRows, fp_csv_over);
            fp_png_over = fullfile(fp_fig_over, pipeline.sanitize_filename(sprintf('peakindex_overview_%s_%s.png', tag, A.name)));
            plot_peakindex_overview(SumRows, metrics, tag, A.name, fp_png_over);
        end
    catch ME
        warning('analyze_peakindex_invertedu: overview failed (%s): %s', A.name, ME.message);
    end
end

end

function plot_peakindex_figure(K, dv, tag, analysisName, fp_png)
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[90 90 900 380]);

levelsCx = categories(K.Complexity);
if isempty(levelsCx)
    levelsCx = {'ComplexityLow','ComplexityHigh'};
end

for i=1:min(2,numel(levelsCx))
    cx = string(levelsCx{i});
    ax = subplot(1,2,i); hold(ax,'on');

    D = K(string(K.Complexity)==cx, :);
    groups = ["Low","High"];
    cols = [0.2 0.5 0.9; 0.9 0.4 0.2];

    for gi=1:2
        g = groups(gi);
        y = D.PeakIndex(string(D.Group)==g);
        if isempty(y); continue; end
        x = gi + (rand(size(y))-0.5)*0.14;
        scatter(ax, x, y, 24, cols(gi,:), 'filled', 'MarkerFaceAlpha',0.6, 'HandleVisibility','off');
        mu = mean(y,'omitnan');
        se = std(y,'omitnan')/sqrt(numel(y));
        errorbar(ax, gi, mu, se, 'o', 'Color', cols(gi,:), 'MarkerFaceColor', cols(gi,:), 'LineWidth',1.4, 'DisplayName', char(g));
    end

    xlim(ax,[0.5 2.5]);
    set(ax,'XTick',1:2,'XTickLabel',{'Low','High'});
    ylabel(ax,'PeakIndex = EEG45 - mean(EEG15, EEG75)');
    title(ax, char(cx), 'Interpreter','none');
    grid(ax,'on');
end

sgtitle(sprintf('Task5 PeakIndex | %s | %s [%s]', analysisName, dv, tag), 'Interpreter','none');
pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function [nsub, mu_all, se_all, p0, bG, pG, bC, pC, bCG, pCG] = extract_peakindex_summary(K, lme)
nsub = numel(unique(string(K.Subject)));
mu_all = mean(K.PeakIndex,'omitnan');
se_all = std(K.PeakIndex,'omitnan')/sqrt(sum(~isnan(K.PeakIndex)));
p0 = NaN; bG = NaN; pG = NaN; bC = NaN; pC = NaN; bCG = NaN; pCG = NaN;

try
    [~, p0] = ttest(K.PeakIndex, 0);
catch
end

try
    C = lme.Coefficients;
    rn = string(C.Name);
    iG = find(contains(rn,'Group'),1,'first');
    iC = find(contains(rn,'Complexity'),1,'first');
    iCG = find(contains(rn,'Complexity') & contains(rn,'Group'),1,'first');
    if ~isempty(iG),  bG = double(C.Estimate(iG)); pG = double(C.pValue(iG)); end
    if ~isempty(iC),  bC = double(C.Estimate(iC)); pC = double(C.pValue(iC)); end
    if ~isempty(iCG), bCG = double(C.Estimate(iCG)); pCG = double(C.pValue(iCG)); end

    try
        A = anova(lme,'DFMethod','Satterthwaite');
        if all(ismember({'Term','pValue'}, A.Properties.VariableNames))
            tt = string(A.Term);
            jG = find(contains(tt,'Group') & ~contains(tt,':'),1,'first');
            jC = find(contains(tt,'Complexity') & ~contains(tt,':'),1,'first');
            jCG = find(contains(tt,'Complexity') & contains(tt,'Group') & contains(tt,':'),1,'first');
            if ~isempty(jG),  pG = double(A.pValue(jG)); end
            if ~isempty(jC),  pC = double(A.pValue(jC)); end
            if ~isempty(jCG), pCG = double(A.pValue(jCG)); end
        end
    catch
    end
catch
end
end

function plot_peakindex_overview(SumRows, metrics, tag, analysisName, fp_png)
% Paper-friendly task5 overview
% row1: mean PeakIndex (+ p vs 0)
% row2: Group effect
% row3: Complexity effect
% row4: Complexity×Group interaction

M = string(metrics(:));
nC = numel(M);
Z = nan(4,nC); P = nan(4,nC); N = nan(1,nC);
for c=1:nC
    i = find(string(SumRows.metric)==M(c),1,'first');
    if isempty(i), continue; end
    Z(1,c) = double(SumRows.mean_peakindex(i));   P(1,c) = double(SumRows.p_vs_zero(i));
    Z(2,c) = double(SumRows.beta_group(i));       P(2,c) = double(SumRows.p_group(i));
    Z(3,c) = double(SumRows.beta_complexity(i));  P(3,c) = double(SumRows.p_complexity(i));
    Z(4,c) = double(SumRows.beta_cxg(i));         P(4,c) = double(SumRows.p_cxg(i));
    N(c)   = double(SumRows.n_subjects(i));
end

set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[90 90 1200 520]);
ax = axes(fig); hold(ax,'on');
imagesc(ax, Z); set(ax,'YDir','normal'); axis(ax,'tight');
colormap(ax, parula(256));
cb = colorbar(ax); cb.Label.String = 'effect / estimate';
mx = max(abs(Z(:)),[],'omitnan'); if isempty(mx)||~isfinite(mx)||mx==0, mx=0.01; end
caxis(ax,[-mx mx]);
set(ax,'XTick',1:nC,'XTickLabel',cellstr(M));
set(ax,'YTick',1:4,'YTickLabel',{'PeakIndex mean','Group','Complexity','Complexity×Group'});
xtickangle(ax,20);
title(ax, sprintf('Task5 PeakIndex overview | %s [%s]', analysisName, tag), 'Interpreter','none');

for r=1:4
    for c=1:nC
        if isnan(Z(r,c)), continue; end
        p = P(r,c); star='';
        if ~isnan(p)
            if p<0.001, star='***'; elseif p<0.01, star='**'; elseif p<0.05, star='*'; end
        end
        if r==1
            txt = sprintf('%.3g%s\np=%.3g\nN=%d', Z(r,c), star, p, round(N(c)));
        else
            txt = sprintf('β=%.3g%s\np=%.3g', Z(r,c), star, p);
        end
        text(ax,c,r,txt,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
    end
end

pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function write_report(fp_md, dv, tag, analysisName, K, lme)
lines = strings(0,1);
lines(end+1) = sprintf('# Task5 PeakIndex (Inverted-U) | %s | %s | %s', analysisName, dv, tag);
lines(end+1) = '';
lines(end+1) = 'Model: `PeakIndex ~ Complexity*Group + (1|Subject)`';
lines(end+1) = '';

% quick means
try
    lines(end+1) = '## Raw means by Complexity × Group';
    lines(end+1) = 'Complexity | Group | n | mean | sem';
    lines(end+1) = '---|---:|---:|---:|---:';
    [G, cx, gr] = findgroups(string(K.Complexity), string(K.Group));
    n = splitapply(@numel, K.PeakIndex, G);
    mu = splitapply(@(x) mean(x,'omitnan'), K.PeakIndex, G);
    se = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), K.PeakIndex, G);
    for i=1:numel(n)
        lines(end+1) = sprintf('%s|%s|%d|%.6g|%.6g', cx(i), gr(i), n(i), mu(i), se(i));
    end
    lines(end+1) = '';
catch
end

lines(end+1) = '## LMM outputs';
lines(end+1) = '- See tables: `peakindex_lmm_fixed_effects_*.csv` and `peakindex_lmm_anova_*.csv`.';
lines(end+1) = '- Interpret Complexity and Group main effects first, then Complexity×Group.';

try
    C = lme.Coefficients;
    rn = string(C.Name);
    ix = find(contains(rn,'Complexity') & contains(rn,'Group'),1);
    if ~isempty(ix)
        lines(end+1) = sprintf('- Complexity×Group p = %.4g', C.pValue(ix));
    end
catch
end

fid = fopen(fp_md,'w');
for i=1:numel(lines)
    fprintf(fid, '%s\n', lines(i));
end
fclose(fid);
end

function w = normalize_wwr(x)
w = string(x);
w = strtrim(w);
for i=1:numel(w)
    tok = regexp(char(w(i)), '(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        w(i) = string(str2double(tok{1}));
    end
end
ok = ismember(w,["15","45","75"]);
w(~ok) = "";
end

function c = normalize_complexity(x)
c = string(x);
c = strtrim(c);
cl = lower(c);
out = repmat("", numel(c), 1);
out(ismember(cl,["low","0","c0","complexitylow"])) = "ComplexityLow";
out(ismember(cl,["high","1","c1","complexityhigh"])) = "ComplexityHigh";
out(c=="ComplexityLow") = "ComplexityLow";
out(c=="ComplexityHigh") = "ComplexityHigh";
isNum = ~isnan(str2double(c));
out(isNum & str2double(c)==0) = "ComplexityLow";
out(isNum & str2double(c)==1) = "ComplexityHigh";
c = out;
end

function g = normalize_high_low(x)
s = string(x);
s = strtrim(s);
sl = lower(s);
g = repmat("", numel(s), 1);
g(ismember(sl,["high","1","高","h"])) = "High";
g(ismember(sl,["low","0","低","l"])) = "Low";
mask = (s=="High" | s=="Low");
g(mask) = s(mask);
end
