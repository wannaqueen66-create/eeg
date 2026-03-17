function out = analyze_individual_checks(AllScene, fp_sum, cfg, tag)
%ANALYZE_INDIVIDUAL_CHECKS Task7: individual-level checks and outlier audit.
%
% Goal:
% - Provide subject-level distribution checks by Group
% - Provide condition-level checks (WWR x Complexity) by Group
% - Block2 O_alpha focused audit for abnormal points
% - Flag potential outlier-driven effects (MAD rule)
% - Summarize significant effects detected from task4/5/6 exported ANOVA tables
%
% Outputs under:
%   <summary>/analysis-2/task7_individual_checks/
%     tables/<tag>/{experience}/...
%     reports/<tag>/{experience}/...
%     figures/<tag>/{experience}/...

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
    struct('name','experience','gcol','ExperienceGroup') ...
};

[fp_root, ~, ~, ~] = pipeline.get_analysis_task_subdirs(fp_sum, 'task7_individual_checks', tag);
if ~exist(fp_root,'dir'); mkdir(fp_root); end

T0 = AllScene;
if ~ismember('subject_id', T0.Properties.VariableNames)
    warning('analyze_individual_checks: missing subject_id. Skipping.');
    return;
end
T0.subject_id = string(T0.subject_id);

% Normalize common factors if present
if ismember('WWR', T0.Properties.VariableNames)
    T0.WWR = normalize_wwr(T0.WWR);
end
if ismember('Complexity', T0.Properties.VariableNames)
    T0.Complexity = normalize_complexity(T0.Complexity);
end

for ai = 1:numel(analyses)
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

    allOut = table();
    allCond = table();

    for mi = 1:numel(metrics)
        dv = metrics(mi);
        if ~ismember(dv, T0.Properties.VariableNames)
            continue;
        end

        % ---------- subject-level mean by group ----------
        dvc = char(dv);
        T = T0(:, {'subject_id', gcolc, dvc});
        T.Properties.VariableNames = {'subject_id','GroupRaw','EEGraw'};
        T.EEG = double(T.EEGraw);
        T.Group = normalize_high_low(T.GroupRaw);
        T = T(~isnan(T.EEG), :);
        T = T(strlength(string(T.Group))>0, :);
        if height(T) < 20
            continue;
        end

        [G, sid, grp] = findgroups(string(T.subject_id), string(T.Group));
        mu = splitapply(@(x) mean(x,'omitnan'), T.EEG, G);
        Sub = table(sid, grp, mu, 'VariableNames', {'subject_id','Group','EEG_mean'});

        [isOut, zmad] = mad_outlier_by_group(Sub.EEG_mean, string(Sub.Group));
        Sub.is_outlier_mad3 = isOut;
        Sub.z_mad = zmad;
        Sub.metric = repmat(dv, height(Sub), 1);
        Sub.analysis = repmat(string(A.name), height(Sub), 1);
        allOut = [allOut; Sub]; %#ok<AGROW>

        writetable(Sub, fullfile(fp_tbl, sprintf('individual_subject_means_%s_%s.csv', dv, tag)));

        try
            fp_png = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('individual_check_%s_%s_%s.png', tag, A.name, dv)));
            plot_individual_figure(Sub, dv, tag, A.name, fp_png);
        catch
        end

        % ---------- condition-level (WWR x Complexity) ----------
        if ismember('WWR', T0.Properties.VariableNames) && ismember('Complexity', T0.Properties.VariableNames)
            Tc = T0(:, {'subject_id', gcolc, dvc, 'WWR', 'Complexity'});
            Tc.Properties.VariableNames = {'subject_id','GroupRaw','EEGraw','WWR','Complexity'};
            Tc.EEG = double(Tc.EEGraw);
            Tc.Group = normalize_high_low(Tc.GroupRaw);
            Tc = Tc(~isnan(Tc.EEG), :);
            Tc = Tc(strlength(string(Tc.Group))>0 & strlength(string(Tc.WWR))>0 & strlength(string(Tc.Complexity))>0, :);
            if height(Tc) >= 24
                [Gc, sidc, grpc, wwrc, cxc] = findgroups(string(Tc.subject_id), string(Tc.Group), string(Tc.WWR), string(Tc.Complexity));
                muc = splitapply(@(x) mean(x,'omitnan'), Tc.EEG, Gc);
                Csub = table(sidc, grpc, wwrc, cxc, muc, 'VariableNames', {'subject_id','Group','WWR','Complexity','EEG_mean'});
                Csub.metric = repmat(dv, height(Csub), 1);
                Csub.analysis = repmat(string(A.name), height(Csub), 1);
                allCond = [allCond; Csub]; %#ok<AGROW>

                writetable(Csub, fullfile(fp_tbl, sprintf('individual_condition_means_%s_%s.csv', dv, tag)));
                try
                    fp_png2 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('individual_condition_box_%s_%s_%s.png', tag, A.name, dv)));
                    plot_condition_box(Csub, dv, tag, A.name, fp_png2);
                catch
                end
            end
        end
    end

    % Write merged tables
    if height(allOut) > 0
        writetable(allOut, fullfile(fp_tbl, sprintf('individual_checks_merged_%s.csv', tag)));
        OutOnly = allOut(logical(allOut.is_outlier_mad3), :);
        writetable(OutOnly, fullfile(fp_tbl, sprintf('individual_outliers_mad3_%s.csv', tag)));

        % Top-N outlier influence on group difference
        Infl = compute_topN_influence(allOut, 3);
        writetable(Infl, fullfile(fp_tbl, sprintf('individual_topN_influence_%s.csv', tag)));
    else
        Infl = table();
    end

    if height(allCond) > 0
        writetable(allCond, fullfile(fp_tbl, sprintf('individual_condition_means_merged_%s.csv', tag)));
    end

    % Block2 O_alpha focused audit
    Oa = table();
    try
        if ismember('O_alpha', T0.Properties.VariableNames) && ismember('block_id', T0.Properties.VariableNames) && ismember('cycle_in_block', T0.Properties.VariableNames)
            Tb = T0(:, {'subject_id', gcolc, 'O_alpha', 'block_id', 'cycle_in_block'});
            Tb.Properties.VariableNames = {'subject_id','GroupRaw','EEG','block_id','cycle_in_block'};
            Tb.Group = normalize_high_low(Tb.GroupRaw);
            Tb = Tb(~isnan(double(Tb.EEG)) & strlength(string(Tb.Group))>0, :);
            Tb = Tb(double(Tb.block_id)==2 & double(Tb.cycle_in_block)>=1 & double(Tb.cycle_in_block)<=6, :);
            if height(Tb) > 0
                Tb.EEG = double(Tb.EEG);
                [isOut2, z2] = mad_outlier_by_cycle_group(Tb.EEG, double(Tb.cycle_in_block), string(Tb.Group));
                Tb.is_outlier_mad3 = isOut2;
                Tb.z_mad = z2;
                Oa = Tb;
                writetable(Oa, fullfile(fp_tbl, sprintf('block2_oalpha_points_%s.csv', tag)));

                fp_png3 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('block2_oalpha_audit_%s_%s.png', tag, A.name)));
                plot_block2_oalpha_audit(Oa, tag, A.name, fp_png3);
            end
        end
    catch
    end

    % Significant effects scan from task4/5/6 tables
    Sig = scan_significant_effects(fp_sum, tag, A.name);
    try
        writetable(Sig, fullfile(fp_tbl, sprintf('significant_effects_detected_%s.csv', tag)));
    catch
    end

    % audit summary + overview + scorecard
    try
        Sum = build_audit_summary(allOut, Oa, Infl);
        if ~isempty(Sum) && height(Sum) > 0
            writetable(Sum, fullfile(fp_tbl, 'task7_audit_summary.csv'));

            fp_png_over = fullfile(fp_fig, 'task7_audit_overview.png');
            plot_audit_overview(Sum, tag, A.name, fp_png_over);

            fp_png_card = fullfile(fp_fig, 'task7_audit_scorecard.png');
            plot_audit_scorecard(Sum, tag, A.name, fp_png_card);
        end
    catch
    end

    % report
    try
        fp_md = fullfile(fp_rep, sprintf('individual_checks_report_%s.md', tag));
        write_report(fp_md, A.name, tag, allOut, allCond, Oa, Sig, Infl);
    catch
    end
end

end

function [isOut, zmad] = mad_outlier_by_group(y, grp)
isOut = false(numel(y),1);
zmad = nan(numel(y),1);
for g = ["Low","High"]
    idx = grp==g;
    yg = y(idx);
    if numel(yg) < 4
        continue;
    end
    med = median(yg,'omitnan');
    madv = mad(yg,1);
    if madv <= eps
        madv = std(yg,'omitnan');
        if madv <= eps; madv = 1; end
    end
    z = abs(yg-med)./madv;
    tmp = false(size(yg));
    tmp(z>3) = true;
    isOut(idx) = tmp;
    zmad(idx) = z;
end
end

function [isOut, zmad] = mad_outlier_by_cycle_group(y, cycle, grp)
isOut = false(numel(y),1);
zmad = nan(numel(y),1);
for c = 1:6
    for g = ["Low","High"]
        idx = (cycle==c) & (grp==g);
        yg = y(idx);
        if numel(yg) < 4
            continue;
        end
        med = median(yg,'omitnan');
        madv = mad(yg,1);
        if madv <= eps
            madv = std(yg,'omitnan');
            if madv <= eps; madv = 1; end
        end
        z = abs(yg-med)./madv;
        tmp = false(size(yg));
        tmp(z>3) = true;
        isOut(idx) = tmp;
        zmad(idx) = z;
    end
end
end

function Infl = compute_topN_influence(allOut, N)
Infl = table('Size',[0 7], ...
    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
    'VariableNames', {'metric','N','diff_full_HighMinusLow','diff_dropTopN','delta','n_total','n_dropped'});

if nargin<2; N=3; end
mets = unique(string(allOut.metric));
for m = mets'
    X = allOut(string(allOut.metric)==m, :);
    if height(X) < 8
        continue;
    end
    muH = mean(X.EEG_mean(string(X.Group)=="High"), 'omitnan');
    muL = mean(X.EEG_mean(string(X.Group)=="Low"), 'omitnan');
    d0 = muH - muL;

    z = abs(double(X.z_mad));
    valid = isfinite(z);
    keep = true(height(X),1);
    k = 0;
    if any(valid)
        idxValid = find(valid);
        [~, ord] = sort(z(valid), 'descend');
        idxDrop = idxValid(ord(1:min(N, numel(ord))));
        k = numel(idxDrop);
        keep(idxDrop) = false;
    end
    Y = X(keep,:);

    muH2 = mean(Y.EEG_mean(string(Y.Group)=="High"), 'omitnan');
    muL2 = mean(Y.EEG_mean(string(Y.Group)=="Low"), 'omitnan');
    d1 = muH2 - muL2;

    Infl = [Infl; {m, N, d0, d1, d1-d0, height(X), k}]; %#ok<AGROW>
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

    yy = Sub.EEG_mean(string(Sub.Group)==g & logical(Sub.is_outlier_mad3));
    if ~isempty(yy)
        xx = gi + zeros(size(yy));
        scatter(ax, xx, yy, 60, 'r', 'o', 'LineWidth',1.2, 'HandleVisibility','off');
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

function plot_condition_box(Csub, dv, tag, analysisName, fp_png)
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[70 70 1200 420]);

cxLevels = ["ComplexityLow","ComplexityHigh"];
wwrLevels = ["15","45","75"];
cols = [0.2 0.5 0.9; 0.9 0.4 0.2];

for ci=1:2
    cx = cxLevels(ci);
    ax = subplot(1,2,ci); hold(ax,'on');
    title(ax, char(cx));

    % explicit legend proxies to avoid auto-binding to wrong box objects
    hLow = scatter(ax, nan, nan, 26, cols(1,:), 'filled', 'MarkerFaceAlpha',0.55, 'DisplayName','Low');
    hHigh = scatter(ax, nan, nan, 26, cols(2,:), 'filled', 'MarkerFaceAlpha',0.55, 'DisplayName','High');

    xbase = 1:3;
    groupOffset = 0.24;
    boxWidth = 0.16;
    jitterWidth = 0.05;
    for gi=1:2
        g = ["Low","High"];
        gg = g(gi);
        for wi=1:3
            w = wwrLevels(wi);
            y = Csub.EEG_mean(string(Csub.Complexity)==cx & string(Csub.WWR)==w & string(Csub.Group)==gg);
            if isempty(y); continue; end
            x0 = xbase(wi) + (gi-1.5)*groupOffset;
            try
                boxchart(ax, repmat(x0, numel(y),1), y, 'BoxWidth',boxWidth, 'MarkerStyle','none', 'BoxFaceColor', cols(gi,:), 'BoxFaceAlpha',0.25, 'HandleVisibility','off');
            catch
                % fallback if boxchart unavailable
                q = quantile(y,[0.25 0.5 0.75]);
                plot(ax,[x0-boxWidth/2 x0+boxWidth/2],[q(2) q(2)],'-','Color',cols(gi,:),'LineWidth',2, 'HandleVisibility','off');
            end
            xj = x0 + (rand(size(y))-0.5)*jitterWidth;
            scatter(ax, xj, y, 14, cols(gi,:), 'filled', 'MarkerFaceAlpha',0.55, 'HandleVisibility','off');
        end
    end

    xlim(ax, [0.5 3.5]);
    set(ax,'XTick',1:3,'XTickLabel',cellstr(wwrLevels));
    xlabel(ax,'WWR');
    ylabel(ax, dv);
    grid(ax,'on');
    legend(ax,[hLow hHigh], 'Location','best');
end

sgtitle(sprintf('Task7 condition-level box/scatter | %s | %s [%s]', analysisName, dv, tag), 'Interpreter','none');
pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function plot_block2_oalpha_audit(Tb, tag, analysisName, fp_png)
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[90 90 980 420]);
ax = axes(fig); hold(ax,'on');

cols = [0.2 0.5 0.9; 0.9 0.4 0.2];
groups = ["Low","High"];

% legend proxy handles
hLow = scatter(ax, nan, nan, 26, cols(1,:), 'filled', 'MarkerFaceAlpha',0.7, 'DisplayName','Low points');
hHigh = scatter(ax, nan, nan, 26, cols(2,:), 'filled', 'MarkerFaceAlpha',0.7, 'DisplayName','High points');
hOut = scatter(ax, nan, nan, 64, 'r', 'o', 'LineWidth',1.1, 'DisplayName','MAD>3 outlier');

for gi=1:2
    g = groups(gi);
    for c=1:6
        y = Tb.EEG(string(Tb.Group)==g & double(Tb.cycle_in_block)==c);
        if isempty(y); continue; end
        x0 = c + (gi-1.5)*0.15;
        x = x0 + (rand(size(y))-0.5)*0.08;
        scatter(ax, x, y, 16, cols(gi,:), 'filled', 'MarkerFaceAlpha',0.55, 'HandleVisibility','off');

        yo = Tb.EEG(string(Tb.Group)==g & double(Tb.cycle_in_block)==c & logical(Tb.is_outlier_mad3));
        if ~isempty(yo)
            xo = x0 + zeros(size(yo));
            scatter(ax, xo, yo, 56, 'r', 'o', 'LineWidth',1.1, 'HandleVisibility','off');
        end

        mu = mean(y,'omitnan');
        se = std(y,'omitnan')/sqrt(numel(y));
        errorbar(ax, x0, mu, se, 'o', 'Color', cols(gi,:), 'MarkerFaceColor', cols(gi,:), 'LineWidth',1.2, 'HandleVisibility','off');
    end
end

set(ax,'XTick',1:6,'XTickLabel',{'B2-1','B2-2','B2-3','B2-4','B2-5','B2-6'});
xlabel(ax,'Block2 cycle');
ylabel(ax,'O_alpha');
grid(ax,'on');
title(ax, sprintf('Task7 Block2 O_alpha audit | %s [%s]', analysisName, tag), 'Interpreter','none');
legend(ax, [hLow hHigh hOut], 'Location','northeastoutside');

pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function Sig = scan_significant_effects(fp_sum, tag, analysisName)
% scan anova tables from task4/5/6 for p<0.05
Sig = table('Size',[0 5], ...
    'VariableTypes', {'string','string','string','string','double'}, ...
    'VariableNames', {'task','table_file','term','p_col','p_value'});

if exist(fullfile(fp_sum,'analysis'),'dir')
    roots = { ...
        fullfile(fp_sum,'analysis','task4_core_lmm_suite',tag,'tables','factor_WWR',analysisName), ...
        fullfile(fp_sum,'analysis','task4_core_lmm_suite',tag,'tables','trend_WWR',analysisName), ...
        fullfile(fp_sum,'analysis','task5_peakindex_invertedu',tag,'tables',analysisName), ...
        fullfile(fp_sum,'analysis','task6_coremetric_special',tag,'tables',analysisName), ...
        fullfile(fp_sum,'analysis','task6_obeta_special',tag,'tables',analysisName) ...
    };
else
    roots = { ...
        fullfile(fp_sum,'analysis-2','task4_core_lmm_suite','factor_WWR','tables',tag,analysisName), ...
        fullfile(fp_sum,'analysis-2','task4_core_lmm_suite','trend_WWR','tables',tag,analysisName), ...
        fullfile(fp_sum,'analysis-2','task5_peakindex_invertedU','tables',tag,analysisName), ...
        fullfile(fp_sum,'analysis-2','task6_obeta_special','tables',tag,analysisName) ...
    };
end

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
elseif contains(s,'task5_peakindex_invertedu')
    tk = "task5";
elseif contains(s,'task6_coremetric_special') || contains(s,'task6_obeta_special')
    tk = "task6";
else
    tk = "unknown";
end
end

function Sum = build_audit_summary(allOut, Oa, Infl)
% Build a compact summary table for paper-style robustness overview.
% Columns: metric, outlier_rate(%), block2_oalpha_outlier_rate(%),
%          influence_abs_pct, direction_flip(flag), robustness_score(A/B/C)

metrics = string(allOut.metric(:));
metrics = unique(metrics,'stable');
if isempty(metrics)
    Sum = table();
    return;
end

m = numel(metrics);
outRate = nan(m,1);
block2Rate = nan(m,1);
inflPct = nan(m,1);
flipFlag = false(m,1);
score = strings(m,1);

% Outlier rate from individual-level check
for i=1:m
    mm = metrics(i);
    idx = string(allOut.metric)==mm;
    if any(idx)
        outRate(i) = 100*mean(double(allOut.is_outlier_mad3(idx)),'omitnan');
    end

    % influence: use first available row for this metric
    ii = find(string(Infl.metric)==mm,1,'first');
    if ~isempty(ii)
        % Preferred: precomputed percentage column (backward/forward compatible)
        if ismember('change_abs_pct', Infl.Properties.VariableNames)
            inflPct(i) = double(Infl.change_abs_pct(ii));
        % Current schema: compute pct from full vs trimmed High-Low differences
        elseif all(ismember({'diff_full_HighMinusLow','diff_dropTopN'}, Infl.Properties.VariableNames))
            d0 = double(Infl.diff_full_HighMinusLow(ii));
            d1 = double(Infl.diff_dropTopN(ii));
            if isfinite(d0) && isfinite(d1)
                inflPct(i) = 100 * abs(d1 - d0) / max(abs(d0), 1e-12);
                flipFlag(i) = (sign(d0)~=sign(d1)) && (d0~=0) && (d1~=0);
            end
        % Legacy schema fallback
        elseif all(ismember({'raw_diff','trimmed_diff'}, Infl.Properties.VariableNames))
            rd = double(Infl.raw_diff(ii));
            td = double(Infl.trimmed_diff(ii));
            if isfinite(rd) && isfinite(td)
                inflPct(i) = 100 * abs(td - rd) / max(abs(rd), 1e-12);
                flipFlag(i) = (sign(rd)~=sign(td)) && (rd~=0) && (td~=0);
            end
        end
    end

    % score rule (simple + interpretable)
    if ~isnan(inflPct(i)) && (inflPct(i) > 40 || flipFlag(i))
        score(i) = "C";
    elseif (~isnan(inflPct(i)) && inflPct(i) > 20) || (~isnan(outRate(i)) && outRate(i) > 12)
        score(i) = "B";
    else
        score(i) = "A";
    end
end

% Block2 O_alpha outlier rate (same value repeated for convenience)
if ~isempty(Oa) && ismember('is_outlier_mad3', Oa.Properties.VariableNames)
    br = 100*mean(double(Oa.is_outlier_mad3),'omitnan');
    block2Rate(:) = br;
end

Sum = table(metrics, outRate, block2Rate, inflPct, flipFlag, score, ...
    'VariableNames', {'metric','outlier_rate_pct','block2_oalpha_outlier_rate_pct','influence_abs_pct','direction_flip','robustness_score'});
end

function plot_audit_overview(Sum, tag, analysisName, fp_png)
% Paper-style overview as a table-like matrix.
% Different metric types use different row-specific color rules; no shared colorbar.

if isempty(Sum) || height(Sum)==0
    return;
end

mets = string(Sum.metric);
n = numel(mets);

outlierVals = double(Sum.outlier_rate_pct(:))';
inflVals = double(Sum.influence_abs_pct(:))';
flipVals = double(Sum.direction_flip(:))';
scoreVals = string(Sum.robustness_score(:))';

rowNames = {'Outlier rate (%)','Influence |d| change (%)','Direction flip (0/1)','Robustness grade'};

set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[60 60 1450 700]);
ax = axes(fig,'Position',[0.18 0.20 0.72 0.68]); hold(ax,'on');
set(ax,'XLim',[0.5 n+0.5],'YLim',[0.5 4.5],'YDir','reverse');

maxOut = max([1, outlierVals]);
maxInfl = max([1, inflVals]);
for c = 1:n
    % row 1: outlier
    t = min(1, outlierVals(c)/maxOut);
    color1 = (1-t)*[0.93 0.96 0.99] + t*[0.29 0.47 0.74];
    rectangle(ax,'Position',[c-0.5,0.5,1,1],'FaceColor',color1,'EdgeColor',[0.85 0.85 0.85]);
    text(ax,c,1,sprintf('%.1f%%', outlierVals(c)),'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Times New Roman','FontSize',11,'FontWeight','bold');

    % row 2: influence
    t = min(1, inflVals(c)/maxInfl);
    if inflVals(c) <= eps
        color2 = [0.96 0.96 0.96];
    else
        color2 = (1-t)*[0.92 0.97 0.94] + t*[0.22 0.63 0.53];
    end
    rectangle(ax,'Position',[c-0.5,1.5,1,1],'FaceColor',color2,'EdgeColor',[0.85 0.85 0.85]);
    text(ax,c,2,sprintf('%.1f%%', inflVals(c)),'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Times New Roman','FontSize',11,'FontWeight','bold');

    % row 3: direction flip
    if flipVals(c) >= 0.5
        color3 = [0.86 0.43 0.38];
        txt3 = '1';
    else
        color3 = [0.94 0.94 0.94];
        txt3 = '0';
    end
    rectangle(ax,'Position',[c-0.5,2.5,1,1],'FaceColor',color3,'EdgeColor',[0.85 0.85 0.85]);
    text(ax,c,3,txt3,'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Times New Roman','FontSize',11,'FontWeight','bold');

    % row 4: robustness
    if scoreVals(c)=="A"
        color4 = [0.79 0.89 0.81];
    elseif scoreVals(c)=="B"
        color4 = [0.97 0.91 0.73];
    else
        color4 = [0.93 0.78 0.78];
    end
    rectangle(ax,'Position',[c-0.5,3.5,1,1],'FaceColor',color4,'EdgeColor',[0.85 0.85 0.85]);
    text(ax,c,4,char(scoreVals(c)),'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Times New Roman','FontSize',11,'FontWeight','bold');
end

set(ax,'XTick',1:n,'XTickLabel',cellstr(mets),'YTick',1:4,'YTickLabel',rowNames);
ax.FontName = 'Times New Roman';
ax.FontSize = 12;
ax.LineWidth = 0.8;
ax.TickLength = [0 0];
box(ax,'on');
grid(ax,'off');
title(ax, sprintf('Task7 audit overview | %s [%s]', analysisName, tag), 'Interpreter','none', 'FontName','Times New Roman', 'FontWeight','normal', 'FontSize',14);

annotation(fig,'textbox',[0.16 0.05 0.70 0.04], ...
    'String','Rule: A = stable; B = moderate sensitivity; C = high sensitivity.', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontName','Times New Roman','FontSize',10);

pipeline.export_figure_png(fig, fp_png, 600);
try; close(fig); catch; end
end

function plot_audit_scorecard(Sum, tag, analysisName, fp_png)
% Paper-style scorecard using panel plots rather than duplicating the overview heatmap.

if isempty(Sum) || height(Sum)==0
    return;
end

mets = string(Sum.metric);
n = numel(mets);
x = 1:n;

sc = nan(n,1);
for i=1:n
    if Sum.robustness_score(i)=="A", sc(i)=3;
    elseif Sum.robustness_score(i)=="B", sc(i)=2;
    elseif Sum.robustness_score(i)=="C", sc(i)=1;
    else, sc(i)=nan;
    end
end

set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[60 60 1500 980]);
tl = tiledlayout(fig,2,2,'TileSpacing','loose','Padding','loose');

ax1 = nexttile(tl,1); hold(ax1,'on');
bar(ax1, x, Sum.outlier_rate_pct, 0.62, 'FaceColor',[0.31 0.47 0.67], 'EdgeColor','none');
ylabel(ax1,'Outlier rate (%)','FontSize',12);
title(ax1,'(a) Outlier rate','FontWeight','normal','FontSize',13);

ax2 = nexttile(tl,2); hold(ax2,'on');
bar(ax2, x, Sum.influence_abs_pct, 0.62, 'FaceColor',[0.20 0.63 0.55], 'EdgeColor','none');
ylabel(ax2,'Influence |d| change (%)','FontSize',12);
title(ax2,'(b) Influence change','FontWeight','normal','FontSize',13);

ax3 = nexttile(tl,3); hold(ax3,'on');
stem(ax3, x, double(Sum.direction_flip), 'filled', 'Color',[0.72 0.34 0.26], 'LineWidth',1.1, 'MarkerSize',6);
ylabel(ax3,'Direction flip','FontSize',12);
ylim(ax3,[-0.05 1.05]);
yticks(ax3,[0 1]);
title(ax3,'(c) Direction flip','FontWeight','normal','FontSize',13);

ax4 = nexttile(tl,4); hold(ax4,'on');
bar(ax4, x, sc, 0.62, 'FaceColor',[0.45 0.45 0.70], 'EdgeColor','none');
ylabel(ax4,'Robustness grade','FontSize',12);
ylim(ax4,[0.5 3.5]);
yticks(ax4,[1 2 3]);
yticklabels(ax4,{'C','B','A'});
title(ax4,'(d) Robustness score','FontWeight','normal','FontSize',13);

allAxes = [ax1 ax2 ax3 ax4];
for k = 1:numel(allAxes)
    ax = allAxes(k);
    set(ax,'XTick',x,'XTickLabel',cellstr(mets));
    xtickangle(ax,0);
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    ax.LineWidth = 0.8;
    box(ax,'off');
    grid(ax,'on');
    ax.GridAlpha = 0.08;
    ax.GridColor = [0 0 0];
end

set([ax1 ax2], 'XTickLabel', []);

for i=1:n
    text(ax1, x(i), Sum.outlier_rate_pct(i) + max(0.3, 0.04*max(Sum.outlier_rate_pct)), sprintf('%.1f%%', Sum.outlier_rate_pct(i)), ...
        'HorizontalAlignment','center','FontSize',10,'FontName','Times New Roman');

    yi = Sum.influence_abs_pct(i);
    textY2 = yi + max(0.12, 0.06*max([1; Sum.influence_abs_pct]));
    if yi <= eps
        textY2 = 0.05;
    end
    text(ax2, x(i), textY2, sprintf('%.1f%%', yi), ...
        'HorizontalAlignment','center','FontSize',10,'FontName','Times New Roman');

    text(ax3, x(i), double(Sum.direction_flip(i)) + 0.05, sprintf('%d', round(double(Sum.direction_flip(i)))), ...
        'HorizontalAlignment','center','FontSize',10,'FontName','Times New Roman');
    text(ax4, x(i), sc(i) + 0.08, char(Sum.robustness_score(i)), ...
        'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','FontName','Times New Roman');
end

xlabel(ax3,'Metric');
xlabel(ax4,'Metric');
title(tl, sprintf('Task7 audit scorecard | %s [%s]', analysisName, tag), 'Interpreter','none', 'FontName','Times New Roman', 'FontWeight','normal', 'FontSize',14);

pipeline.export_figure_png(fig, fp_png, 600);
try; close(fig); catch; end
end

function write_report(fp_md, analysisName, tag, allOut, allCond, Oa, Sig, Infl)
lines = strings(0,1);
lines(end+1) = sprintf('# Task7 Individual checks | %s | %s', analysisName, tag);
lines(end+1) = '';
lines(end+1) = 'This report checks whether significant effects may be driven by a few extreme subjects.';
lines(end+1) = '';

if ~isempty(allOut) && height(allOut)>0
    nAll = height(allOut);
    nOut = sum(logical(allOut.is_outlier_mad3));
    lines(end+1) = sprintf('- Total subject-level points: %d', nAll);
    lines(end+1) = sprintf('- MAD>3 outliers flagged: %d', nOut);
end
if ~isempty(allCond) && height(allCond)>0
    lines(end+1) = sprintf('- Condition-level subject means (WWR×Complexity) rows: %d', height(allCond));
end
if ~isempty(Oa) && height(Oa)>0
    lines(end+1) = sprintf('- Block2 O_alpha points audited: %d', height(Oa));
    lines(end+1) = sprintf('- Block2 O_alpha MAD>3 points: %d', sum(logical(Oa.is_outlier_mad3)));
end
lines(end+1) = '';

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

if ~isempty(Infl) && height(Infl)>0
    lines(end+1) = '';
    lines(end+1) = '## Top-N outlier influence on High-Low difference';
    lines(end+1) = 'metric | diff_full | diff_dropTopN | delta';
    lines(end+1) = '---|---:|---:|---:';
    for i=1:height(Infl)
        lines(end+1) = sprintf('%s|%.6g|%.6g|%.6g', string(Infl.metric(i)), Infl.diff_full_HighMinusLow(i), Infl.diff_dropTopN(i), Infl.delta(i));
    end
end

lines(end+1) = '';
lines(end+1) = '## Files';
lines(end+1) = '- tables: subject means, condition means, outlier list, significant effects table, topN influence';
lines(end+1) = '- figures: per-metric subject distribution, condition-level box/scatter, Block2 O_alpha audit';

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
