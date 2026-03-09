function out = analyze_obeta_special(AllScene, fp_sum, cfg, tag)
%ANALYZE_OBETA_SPECIAL Task6: core-metric special robustness models.
%
% Default metrics are aligned with the repository's main core metric set:
%   O_alpha, O_theta, O_beta, F_theta
%
% For each grouping line (ExperienceGroup / SportFreqGroup), fit:
%   Model A: EEG ~ Group + (1|Subject)
%   Model B: EEG ~ WWR + Complexity + Group + (1|Subject)
%
% Outputs under:
%   task6_obeta_special (historical task name retained for compatibility)
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
    warning('analyze_obeta_special: fitlme not found (Stats toolbox missing). Skipping.');
    return;
end

req = {'subject_id'};
for i=1:numel(req)
    if ~ismember(req{i}, AllScene.Properties.VariableNames)
        warning('analyze_obeta_special: missing required column %s. Skipping.', req{i});
        return;
    end
end

metrics = ["O_alpha","O_theta","O_beta","F_theta"];
try
    if isfield(cfg,'task6_metrics') && ~isempty(cfg.task6_metrics)
        metrics = string(cfg.task6_metrics);
    end
catch
end

hasAnyMetric = any(ismember(metrics, string(AllScene.Properties.VariableNames)));
if ~hasAnyMetric
    warning('analyze_obeta_special: none of task6_metrics are present. Skipping.');
    return;
end

analyses = { ...
    struct('name','experience','gcol','ExperienceGroup'), ...
    struct('name','sportfreq','gcol','SportFreqGroup') ...
};

[fp_root, ~, ~, ~] = pipeline.get_analysis_task_subdirs(fp_sum, 'task6_obeta_special', tag);

T0 = AllScene;
T0.subject_id = string(T0.subject_id);
if ismember('WWR', T0.Properties.VariableNames)
    T0.WWR = normalize_wwr(T0.WWR);
end
if ismember('Complexity', T0.Properties.VariableNames)
    T0.Complexity = normalize_complexity(T0.Complexity);
end

for ai=1:numel(analyses)
    A = analyses{ai};
    gcol = string(A.gcol);
    if ~ismember(gcol, string(T0.Properties.VariableNames))
        continue;
    end
    gcolc = char(gcol);

    fp_tbl = fullfile(fp_root, 'tables', tag, A.name);
    fp_rep = fullfile(fp_root, 'reports', tag, A.name);
    fp_fig = fullfile(fp_root, 'figures', tag, A.name);
    if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
    if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
    if ~exist(fp_fig,'dir'); mkdir(fp_fig); end

    for mi=1:numel(metrics)
        met = string(metrics(mi));
        if ~ismember(met, string(T0.Properties.VariableNames))
            continue;
        end
        metc = char(met);

        vars = {'subject_id',metc,gcolc};
        if ismember('WWR', T0.Properties.VariableNames)
            vars{end+1} = 'WWR';
        end
        if ismember('Complexity', T0.Properties.VariableNames)
            vars{end+1} = 'Complexity';
        end

        T = T0(:, vars);
        T.Properties.VariableNames{3} = 'GroupRaw';
        T.EEG = double(T.(metc));
        T.Group = normalize_high_low(T.GroupRaw);
        T = T(~isnan(T.EEG), :);
        T = T(strlength(string(T.Group))>0, :);
        if height(T) < 20
            continue;
        end

        T.Subject = categorical(string(T.subject_id));
        T.Group = categorical(string(T.Group), {'Low','High'});

        hasControls = ismember('WWR', T.Properties.VariableNames) && ismember('Complexity', T.Properties.VariableNames);
        if ismember('WWR', T.Properties.VariableNames)
            T.WWR = categorical(string(T.WWR), {'15','45','75'});
        end
        if ismember('Complexity', T.Properties.VariableNames)
            T.Complexity = categorical(string(T.Complexity), {'ComplexityLow','ComplexityHigh'});
        end

        % For fair Model A vs Model B comparison, use the same complete-case sample
        % when control predictors (WWR, Complexity) are available.
        if hasControls
            useCtrl = ~isundefined(T.WWR) & ~isundefined(T.Complexity);
            T = T(useCtrl, :);
            if height(T) < 20
                continue;
            end
        end

        % Save analysis-ready
        keep = {'Subject','Group','EEG'};
        if ismember('WWR', T.Properties.VariableNames); keep{end+1} = 'WWR'; end
        if ismember('Complexity', T.Properties.VariableNames); keep{end+1} = 'Complexity'; end
        writetable(T(:,keep), fullfile(fp_tbl, sprintf('%s_analysis_ready_%s.csv', lower(metc), tag)));

        % Model A
        try
            lmeA = fitlme(T, 'EEG ~ Group + (1|Subject)');
        catch ME
            warning('analyze_obeta_special: Model A failed (%s/%s/%s): %s', tag, A.name, met, ME.message);
            continue;
        end

        write_lme_tables(lmeA, fp_tbl, sprintf('%s_modelA_group_only_%s', lower(metc), tag));

        % Model B (controlled)
        lmeB = [];
        try
            if ismember('WWR', T.Properties.VariableNames) && ismember('Complexity', T.Properties.VariableNames)
                lmeB = fitlme(T, 'EEG ~ WWR + Complexity + Group + (1|Subject)');
                write_lme_tables(lmeB, fp_tbl, sprintf('%s_modelB_controlled_%s', lower(metc), tag));
            end
        catch ME
            warning('analyze_obeta_special: Model B failed (%s/%s/%s): %s', tag, A.name, met, ME.message);
        end

        % Figure 1: raw group scatter/mean
        try
            fp_png = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('%s_special_%s_%s.png', lower(metc), tag, A.name)));
            plot_group_metric(T, metc, tag, A.name, fp_png);
        catch
        end

        % Figure 2: Model A vs Model B group-effect robustness (Estimate + p)
        try
            fp_png2 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('%s_group_model_compare_%s_%s.png', lower(metc), tag, A.name)));
            plot_group_model_compare(lmeA, lmeB, metc, tag, A.name, fp_png2);
        catch
        end

        % Report
        try
            fp_md = fullfile(fp_rep, sprintf('%s_special_report_%s.md', lower(metc), tag));
            write_report(fp_md, metc, tag, A.name, T, lmeA, lmeB);
        catch
        end
    end
end

end

function write_lme_tables(lme, fp_tbl, stem)
try
    C = to_table_compat(lme.Coefficients);
    try
        C = pipeline.add_holm_to_anova_table(C);
    catch
    end
    writetable(C, fullfile(fp_tbl, sprintf('%s_fixed_effects.csv', stem)));
catch ME
    warning('analyze_obeta_special: failed to write fixed effects (%s): %s', stem, ME.message);
end
try
    A = to_table_compat(anova(lme,'DFMethod','Satterthwaite'));
    try
        A = pipeline.add_holm_to_anova_table(A);
    catch
    end
    writetable(A, fullfile(fp_tbl, sprintf('%s_anova.csv', stem)));
catch ME
    warning('analyze_obeta_special: failed to write ANOVA (%s): %s', stem, ME.message);
end
end

function plot_group_metric(T, metricName, tag, analysisName, fp_png)
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 760 360]);
ax = axes(fig); hold(ax,'on');

groups = ["Low","High"];
cols = [0.2 0.5 0.9; 0.9 0.4 0.2];
for gi=1:2
    g = groups(gi);
    y = T.EEG(string(T.Group)==g);
    if isempty(y); continue; end
    x = gi + (rand(size(y))-0.5)*0.14;
    scatter(ax, x, y, 22, cols(gi,:), 'filled', 'MarkerFaceAlpha',0.6, 'HandleVisibility','off');
    mu = mean(y,'omitnan');
    se = std(y,'omitnan')/sqrt(numel(y));
    errorbar(ax, gi, mu, se, 'o', 'Color', cols(gi,:), 'MarkerFaceColor', cols(gi,:), 'LineWidth',1.4, 'DisplayName', char(g));
end
set(ax,'XTick',1:2,'XTickLabel',{'Low','High'});
xlim(ax,[0.5 2.5]);
ylabel(ax, strrep(metricName,'_','\\_'), 'Interpreter','none');
grid(ax,'on');
title(ax, sprintf('Task6 %s special | %s [%s]', strrep(metricName,'_','\\_'), analysisName, tag), 'Interpreter','none');
legend(ax,'Location','best');
pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function plot_group_model_compare(lmeA, lmeB, metricName, tag, analysisName, fp_png)
% Visual summary for task6 robustness check:
% - left: Group(High vs Low) estimate ±95%CI in Model A/B
% - right: Group p-values in Model A/B as -log10(p), with p=0.05 line

if isempty(lmeA)
    return;
end

% Extract Group term from fixed-effects tables
[eA, seA, pA] = extract_group_term(lmeA);
if isnan(eA)
    return;
end

hasB = ~isempty(lmeB);
if hasB
    [eB, seB, pB] = extract_group_term(lmeB);
else
    eB = NaN; seB = NaN; pB = NaN;
end

set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 920 360]);

% Left panel: estimate + 95% CI
ax1 = subplot(1,2,1); hold(ax1,'on');
ciA = [eA - 1.96*seA, eA + 1.96*seA];
errorbar(ax1, 1, eA, eA-ciA(1), ciA(2)-eA, 'o', 'LineWidth',1.6, ...
    'MarkerFaceColor',[0.2 0.5 0.9], 'Color',[0.2 0.5 0.9]);
if hasB && isfinite(eB)
    ciB = [eB - 1.96*seB, eB + 1.96*seB];
    errorbar(ax1, 2, eB, eB-ciB(1), ciB(2)-eB, 'o', 'LineWidth',1.6, ...
        'MarkerFaceColor',[0.9 0.4 0.2], 'Color',[0.9 0.4 0.2]);
    set(ax1,'XTick',[1 2],'XTickLabel',{'Model A','Model B'});
    xlim(ax1,[0.5 2.5]);
else
    ciB = [NaN NaN];
    set(ax1,'XTick',1,'XTickLabel',{'Model A'});
    xlim(ax1,[0.5 1.5]);
end
yline(ax1, 0, 'k--');
ylabel(ax1, 'Group Estimate (High vs Low)');
title(ax1, sprintf('%s group effect size', strrep(metricName,'_','\\_')), 'Interpreter','none');
grid(ax1,'on');

% Add numeric labels on left panel (beta + 95%CI + p)
text(ax1, 1, eA, sprintf('  β=%.4g\n  95%%CI=[%.4g, %.4g]\n  p=%.4g', eA, ciA(1), ciA(2), pA), ...
    'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8,'Color',[0.1 0.3 0.6]);
if hasB && isfinite(eB)
    text(ax1, 2, eB, sprintf('  β=%.4g\n  95%%CI=[%.4g, %.4g]\n  p=%.4g', eB, ciB(1), ciB(2), pB), ...
        'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8,'Color',[0.6 0.2 0.1]);
end

% Right panel: p-values as -log10(p)
ax2 = subplot(1,2,2); hold(ax2,'on');
if hasB && isfinite(pB)
    vals = -log10([max(pA,realmin), max(pB,realmin)]);
    bar(ax2, [1 2], vals, 0.5);
    set(ax2,'XTick',[1 2],'XTickLabel',{'Model A','Model B'});
    xlim(ax2,[0.5 2.5]);
    text(ax2,1,vals(1),sprintf('  p=%.4g',pA),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8);
    text(ax2,2,vals(2),sprintf('  p=%.4g',pB),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8);
else
    vals = -log10(max(pA,realmin));
    bar(ax2, 1, vals, 0.5);
    set(ax2,'XTick',1,'XTickLabel',{'Model A'});
    xlim(ax2,[0.5 1.5]);
    text(ax2,1,vals,sprintf('  p=%.4g',pA),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8);
end
yline(ax2, -log10(0.05), 'r--', 'p=0.05', 'LabelVerticalAlignment','bottom');
ylabel(ax2, '-log10(p)');
title(ax2, sprintf('Group p-values (A=%.3g%s)', pA, ternary(hasB && isfinite(pB), sprintf(', B=%.3g', pB), '')));
grid(ax2,'on');

sgtitle(sprintf('Task6 %s group robustness | %s [%s]', strrep(metricName,'_','\\_'), analysisName, tag), 'Interpreter','none');
pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function [est, se, p] = extract_group_term(lme)
est = NaN; se = NaN; p = NaN;
try
    C = to_table_compat(lme.Coefficients);
    names = string(C.Name);
    idx = find(contains(lower(names),'group'),1,'first');
    if isempty(idx)
        return;
    end
    est = double(C.Estimate(idx));
    se  = double(C.SE(idx));
    p   = double(C.pValue(idx));

    % Prefer ANOVA Group p when available
    try
        A = to_table_compat(anova(lme,'DFMethod','Satterthwaite'));
        if ismember('Term', A.Properties.VariableNames) && ismember('pValue', A.Properties.VariableNames)
            t = string(A.Term);
            ia = find(contains(lower(t),'group'),1,'first');
            if ~isempty(ia)
                p = double(A.pValue(ia));
            end
        end
    catch
    end
catch
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function write_report(fp_md, metricName, tag, analysisName, T, lmeA, lmeB)
lines = strings(0,1);
lines(end+1) = sprintf('# Task6 %s special | %s | %s', strrep(metricName,'_','\\_'), analysisName, tag);
lines(end+1) = '';
lines(end+1) = 'Models:';
lines(end+1) = '- Model A: `EEG ~ Group + (1|Subject)`';
if ~isempty(lmeB)
    lines(end+1) = '- Model B: `EEG ~ WWR + Complexity + Group + (1|Subject)`';
end
lines(end+1) = '';

% Raw means
try
    [G, grp] = findgroups(string(T.Group));
    n = splitapply(@numel, T.EEG, G);
    mu = splitapply(@(x) mean(x,'omitnan'), T.EEG, G);
    se = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T.EEG, G);
    lines(end+1) = '## Raw group means';
    lines(end+1) = 'Group | n | mean | sem';
    lines(end+1) = '---|---:|---:|---:';
    for i=1:numel(n)
        lines(end+1) = sprintf('%s|%d|%.6g|%.6g', grp(i), n(i), mu(i), se(i));
    end
    lines(end+1) = '';
catch
end

lines(end+1) = '## LMM outputs';
lines(end+1) = '- See exported fixed_effects/anova CSVs in tables/.';
lines(end+1) = '- Compare Group effect in Model A vs Model B to assess robustness after controlling WWR & Complexity.';
lines(end+1) = sprintf('- Visual summary PNG: `%s_group_model_compare_<tag>_<analysis>.png` (estimate±95%%CI + p-value comparison).', lower(metricName));

fid = fopen(fp_md,'w');
for i=1:numel(lines)
    fprintf(fid, '%s\n', lines(i));
end
fclose(fid);
end

function Ttbl = to_table_compat(X)
if istable(X)
    Ttbl = X;
    return;
end
try
    if isa(X,'dataset')
        Ttbl = dataset2table(X);
        return;
    end
catch
end
try
    Ttbl = struct2table(X);
    return;
catch
end
error('Unsupported output type for table export: %s', class(X));
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
