function out = analyze_individual_checks(AllScene, fp_sum, cfg, tag)
%ANALYZE_INDIVIDUAL_CHECKS Task7: individual-level checks and outlier audit.
%
% Goal:
% - Provide subject-level distribution checks (box/scatter) by Group
% - Flag potential outlier-driven effects (MAD rule)
% - Summarize significant effects detected from task4/5/6 exported ANOVA tables
%
% Outputs under:
%   <summary>/analysis-2/task7_individual_checks/
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

metrics = ["O_theta","F_theta","O_alpha","O_beta"];
try
    if isfield(cfg,'task7_metrics') && ~isempty(cfg.task7_metrics)
        metrics = string(cfg.task7_metrics);
    end
catch
end

analyses = { ...
    struct('name','experience','gcol','ExperienceGroup'), ...
    struct('name','sportfreq','gcol','SportFreqGroup') ...
};

fp_root = fullfile(fp_sum, 'analysis-2', 'task7_individual_checks');
if ~exist(fp_root,'dir'); mkdir(fp_root); end

T0 = AllScene;
if ~ismember('subject_id', T0.Properties.VariableNames)
    warning('analyze_individual_checks: missing subject_id. Skipping.');
    return;
end
T0.subject_id = string(T0.subject_id);

for ai = 1:numel(analyses)
    A = analyses{ai};
    gcol = string(A.gcol);
    if ~ismember(gcol, T0.Properties.VariableNames)
        continue;
    end

    fp_tbl = fullfile(fp_root, 'tables', tag, A.name);
    fp_rep = fullfile(fp_root, 'reports', tag, A.name);
    fp_fig = fullfile(fp_root, 'figures', tag, A.name);
    if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
    if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
    if ~exist(fp_fig,'dir'); mkdir(fp_fig); end

    allOut = table();

    for mi = 1:numel(metrics)
        dv = metrics(mi);
        if ~ismember(dv, T0.Properties.VariableNames)
            continue;
        end

        T = T0(:, {'subject_id', gcol, dv});
        T.Properties.VariableNames = {'subject_id','GroupRaw','EEGraw'};
        T.EEG = double(T.EEGraw);
        T.Group = normalize_high_low(T.GroupRaw);
        T = T(~isnan(T.EEG), :);
        T = T(strlength(string(T.Group))>0, :);
        if height(T) < 20
            continue;
        end

        % Subject-level means (one value per subject)
        [G, sid, grp] = findgroups(string(T.subject_id), string(T.Group));
        mu = splitapply(@(x) mean(x,'omitnan'), T.EEG, G);
        Sub = table(sid, grp, mu, 'VariableNames', {'subject_id','Group','EEG_mean'});

        % Outlier check by group: |x - median| > 3*MAD
        isOut = false(height(Sub),1);
        zmad = nan(height(Sub),1);
        for g = ["Low","High"]
            idx = string(Sub.Group)==g;
            y = Sub.EEG_mean(idx);
            if numel(y) < 4
                continue;
            end
            med = median(y,'omitnan');
            madv = mad(y,1); % median absolute deviation (scaled=1)
            if madv <= eps
                madv = std(y,'omitnan');
                if madv <= eps; madv = 1; end
            end
            z = abs(y - med) ./ madv;
            tmp = false(size(y));
            tmp(z > 3) = true;
            isOut(idx) = tmp;
            zmad(idx) = z;
        end
        Sub.is_outlier_mad3 = isOut;
        Sub.z_mad = zmad;
        Sub.metric = repmat(dv, height(Sub), 1);
        Sub.analysis = repmat(string(A.name), height(Sub), 1);
        allOut = [allOut; Sub]; %#ok<AGROW>

        % write per-metric subject table
        writetable(Sub, fullfile(fp_tbl, sprintf('individual_subject_means_%s_%s.csv', dv, tag)));

        % figure
        try
            fp_png = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('individual_check_%s_%s_%s.png', tag, A.name, dv)));
            plot_individual_figure(Sub, dv, tag, A.name, fp_png);
        catch
        end
    end

    % Write merged outlier table
    if height(allOut) > 0
        writetable(allOut, fullfile(fp_tbl, sprintf('individual_checks_merged_%s.csv', tag)));
        OutOnly = allOut(logical(allOut.is_outlier_mad3), :);
        writetable(OutOnly, fullfile(fp_tbl, sprintf('individual_outliers_mad3_%s.csv', tag)));
    end

    % Significant effects scan from task4/5/6 tables
    Sig = scan_significant_effects(fp_sum, tag, A.name);
    try
        writetable(Sig, fullfile(fp_tbl, sprintf('significant_effects_detected_%s.csv', tag)));
    catch
    end

    % report
    try
        fp_md = fullfile(fp_rep, sprintf('individual_checks_report_%s.md', tag));
        write_report(fp_md, A.name, tag, allOut, Sig);
    catch
    end
end

end

function plot_individual_figure(Sub, dv, tag, analysisName, fp_png)
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[100 100 780 380]);
ax = axes(fig); hold(ax,'on');

groups = ["Low","High"];
cols = [0.2 0.5 0.9; 0.9 0.4 0.2];

for gi=1:2
    g = groups(gi);
    y = Sub.EEG_mean(string(Sub.Group)==g);
    if isempty(y); continue; end
    x = gi + (rand(size(y))-0.5)*0.15;
    scatter(ax, x, y, 22, cols(gi,:), 'filled', 'MarkerFaceAlpha',0.6, 'HandleVisibility','off');
    mu = mean(y,'omitnan');
    se = std(y,'omitnan')/sqrt(numel(y));
    errorbar(ax, gi, mu, se, 'o', 'Color', cols(gi,:), 'MarkerFaceColor', cols(gi,:), 'LineWidth',1.4, 'DisplayName', char(g));

    % highlight outliers (red ring)
    yy = Sub.EEG_mean(string(Sub.Group)==g & logical(Sub.is_outlier_mad3));
    if ~isempty(yy)
        xx = gi + zeros(size(yy));
        scatter(ax, xx, yy, 60, 'r', 'o', 'LineWidth',1.2, 'DisplayName', sprintf('%s outlier', g));
    end
end

set(ax,'XTick',1:2,'XTickLabel',{'Low','High'});
xlim(ax,[0.5 2.5]);
ylabel(ax, dv);
grid(ax,'on');
title(ax, sprintf('Task7 individual check | %s | %s [%s]', analysisName, dv, tag), 'Interpreter','none');
legend(ax,'Location','best');

pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function Sig = scan_significant_effects(fp_sum, tag, analysisName)
% scan anova tables from task4/5/6 for p<0.05
Sig = table('Size',[0 5], ...
    'VariableTypes', {'string','string','string','string','double'}, ...
    'VariableNames', {'task','table_file','term','p_col','p_value'});

roots = { ...
    fullfile(fp_sum,'analysis-2','task4_core_lmm_suite','factor_WWR','tables',tag,analysisName), ...
    fullfile(fp_sum,'analysis-2','task4_core_lmm_suite','trend_WWR','tables',tag,analysisName), ...
    fullfile(fp_sum,'analysis-2','task5_peakindex_invertedU','tables',tag,analysisName), ...
    fullfile(fp_sum,'analysis-2','task6_obeta_special','tables',tag,analysisName) ...
};

for ri=1:numel(roots)
    r = roots{ri};
    if ~exist(r,'dir'); continue; end
    D = dir(fullfile(r, '*anova*.csv'));
    for i=1:numel(D)
        fp = fullfile(D(i).folder, D(i).name);
        try
            T = readtable(fp, 'TextType','string');
        catch
            continue;
        end
        if isempty(T) || height(T)==0; continue; end

        [pcol, termcol] = detect_cols(T);
        if pcol==""; continue; end
        p = double(T.(pcol));
        for k=1:height(T)
            if ~isnan(p(k)) && p(k) < 0.05
                term = "";
                if termcol ~= ""
                    term = string(T.(termcol)(k));
                end
                tk = infer_task_from_path(fp);
                Sig = [Sig; {tk, string(fp), term, pcol, p(k)}]; %#ok<AGROW>
            end
        end
    end
end
end

function [pcol, termcol] = detect_cols(T)
pcol = ""; termcol = "";
vars = string(T.Properties.VariableNames);
for n = ["pValue","pvalue","p","ProbF","PrF"]
    if any(vars==n)
        pcol = n; break;
    end
end
for n = ["Term","Name","Effect","Source"]
    if any(vars==n)
        termcol = n; break;
    end
end
end

function tk = infer_task_from_path(fp)
s = lower(string(fp));
if contains(s,'task4_core_lmm_suite')
    tk = "task4";
elseif contains(s,'task5_peakindex_invertedu') || contains(s,'task5_peakindex_invertedu')
    tk = "task5";
elseif contains(s,'task6_obeta_special')
    tk = "task6";
else
    tk = "unknown";
end
end

function write_report(fp_md, analysisName, tag, allOut, Sig)
lines = strings(0,1);
lines(end+1) = sprintf('# Task7 Individual checks | %s | %s', analysisName, tag);
lines(end+1) = '';
lines(end+1) = 'This report checks whether significant effects may be driven by a few extreme subjects.';
lines(end+1) = '';

if ~isempty(allOut)
    nAll = height(allOut);
    nOut = sum(logical(allOut.is_outlier_mad3));
    lines(end+1) = sprintf('- Total subject-level points: %d', nAll);
    lines(end+1) = sprintf('- MAD>3 outliers flagged: %d', nOut);
    lines(end+1) = '';
end

if isempty(Sig) || height(Sig)==0
    lines(end+1) = '## Significant effects scan';
    lines(end+1) = '- No p<0.05 terms detected from task4/task5/task6 ANOVA files (or files not found).';
else
    lines(end+1) = '## Significant effects scan (p<0.05)';
    lines(end+1) = 'task | term | p';
    lines(end+1) = '---|---|---:';
    for i=1:height(Sig)
        lines(end+1) = sprintf('%s|%s|%.4g', string(Sig.task(i)), string(Sig.term(i)), Sig.p_value(i));
    end
end

lines(end+1) = '';
lines(end+1) = '## Files';
lines(end+1) = '- tables: subject means, merged checks, outlier list, significant effects table';
lines(end+1) = '- figures: per-metric distribution plots with outlier markers';

fid = fopen(fp_md,'w');
for i=1:numel(lines)
    fprintf(fid, '%s\n', lines(i));
end
fclose(fid);
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
