function fp = plot_trialindex_trend_facets(Tready, analysisName, metric, tag, fp_out, lme)
%PLOT_TRIALINDEX_TREND_FACETS Faceted paper plot by WWR x Complexity.
%
% Creates a 2 (ComplexityLow/High) x 3 (WWR 15/45/75) grid.
% Within each panel: mean±SEM of EEG by scene-specific occurrence index
% (1st/2nd/...) for Group Low/High.
%
% Why occurrence index (not global TrialIndex 1..12)?
%   In a single scene panel (e.g., WWR45), global TrialIndex mixes different
%   subject subsets at each x and can confound sequence/order assignment.
%   Occurrence index aligns "same scene, kth exposure" within subject.
%
% Requires Tready columns: subject_id, Group, WWR, Complexity, TrialIndex, EEG

if nargin < 3; metric = "EEG"; end
if nargin < 4; tag = "raw"; end
if nargin < 5 || strlength(string(fp_out))==0
    fp_out = pwd;
end

if ~exist(fp_out,'dir'); mkdir(fp_out); end

set(0,'DefaultFigureVisible','off');

% normalize to string
grp = string(Tready.Group);
wwr = string(Tready.WWR);
cx  = string(Tready.Complexity);
xTrial = double(Tready.TrialIndex);
y = double(Tready.EEG);

if ismember('subject_id', Tready.Properties.VariableNames)
    sid = string(Tready.subject_id);
elseif ismember('Subject', Tready.Properties.VariableNames)
    sid = string(Tready.Subject);
else
    sid = repmat("", height(Tready), 1);
end

% keep High/Low and non-missing
mask = (grp=="High" | grp=="Low") & ~isnan(xTrial) & ~isnan(y) & ...
       strlength(wwr)>0 & strlength(cx)>0 & strlength(sid)>0;
grp = grp(mask); wwr = wwr(mask); cx = cx(mask); xTrial = xTrial(mask); y = y(mask); sid = sid(mask);

if numel(y) < 40
    fp = "";
    return;
end

% preferred orders
wwrLevels = ["15","45","75"];
cxLevels  = ["ComplexityLow","ComplexityHigh"];

% if WWR stored as numbers but as strings already ok; keep only those
okW = ismember(wwr, wwrLevels);
okC = ismember(cx, cxLevels);
if ~any(okW)
    % fallback: try extract digits
    w2 = repmat("", numel(wwr), 1);
    for i=1:numel(wwr)
        tok = regexp(char(wwr(i)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            w2(i) = string(str2double(tok{1}));
        end
    end
    wwr = w2;
    okW = ismember(wwr, wwrLevels);
end
if ~any(okC)
    % fallback: accept low/high
    c2 = lower(strtrim(cx));
    c3 = repmat("", numel(c2), 1);
    c3(c2=="low" | c2=="c0") = "ComplexityLow";
    c3(c2=="high"| c2=="c1") = "ComplexityHigh";
    cx = c3;
    okC = ismember(cx, cxLevels);
end

keep = okW & okC;
grp = grp(keep); wwr = wwr(keep); cx = cx(keep); xTrial = xTrial(keep); y = y(keep); sid = sid(keep);
if numel(y) < 40
    fp = "";
    return;
end

colors = struct('Low',[0.2 0.5 0.9], 'High',[0.9 0.4 0.2]);

fig = figure('Position',[80 80 1180 680], 'Color','w');

for r=1:numel(cxLevels)
    for c=1:numel(wwrLevels)
        cxv = cxLevels(r);
        wv  = wwrLevels(c);
        ax = subplot(numel(cxLevels), numel(wwrLevels), (r-1)*numel(wwrLevels)+c);
        hold(ax,'on');

        idxPanel = (cx==cxv & wwr==wv);
        if sum(idxPanel) < 10
            title(ax, sprintf('%s | WWR%s (n<10)', cxv, wv), 'Interpreter','none');
            grid(ax,'on');
            continue;
        end

        % Build panel-local scene occurrence index per subject:
        % sort by global TrialIndex and assign 1..k within subject.
        panelRows = find(idxPanel);
        sidP = sid(panelRows);
        trP  = xTrial(panelRows);

        [~, ord] = sortrows([double(categorical(sidP)) trP], [1 2]);
        panelRowsSorted = panelRows(ord);
        sidSorted = sid(panelRowsSorted);

        occSorted = nan(numel(panelRowsSorted),1);
        u = unique(sidSorted);
        for ui=1:numel(u)
            ii = find(sidSorted==u(ui));
            occSorted(ii) = (1:numel(ii))';
        end

        % map back to panel row order
        occ = nan(numel(panelRows),1);
        occ(ord) = occSorted;

        % hard cap in case of malformed duplicates
        occ(~isfinite(occ) | occ<1) = NaN;

        maxOcc = max(occ,[],'omitnan');
        if isempty(maxOcc) || ~isfinite(maxOcc) || maxOcc < 1
            title(ax, sprintf('%s | WWR%s (no valid occ idx)', cxv, wv), 'Interpreter','none');
            grid(ax,'on');
            continue;
        end
        occLevels = (1:floor(maxOcc))';

        for gname = ["Low","High"]
            idxG = idxPanel & (grp==gname);
            if sum(idxG) < 5
                continue;
            end

            rowsG = find(idxG);
            [lia, loc] = ismember(rowsG, panelRows);
            loc = loc(lia);
            rowsG = rowsG(lia);

            occG = occ(loc);
            yG = y(rowsG);

            mu = nan(numel(occLevels),1);
            sem = nan(numel(occLevels),1);
            for i=1:numel(occLevels)
                oi = occLevels(i);
                yy = yG(occG==oi);
                if isempty(yy)
                    continue;
                end
                mu(i) = mean(yy,'omitnan');
                sem(i) = std(yy,'omitnan')/sqrt(numel(yy));
            end

            cc = colors.(char(gname));
            errorbar(ax, occLevels, mu, sem, 'o-', 'LineWidth', 1.4, 'Color', cc, ...
                'MarkerFaceColor', cc, 'DisplayName', char(gname));

            % visual linear fit
            try
                ok = ~isnan(mu);
                if sum(ok) >= 2
                    p = polyfit(occLevels(ok), mu(ok), 1);
                    yhat = polyval(p, occLevels);
                    plot(ax, occLevels, yhat, '--', 'LineWidth', 1.1, 'Color', cc, 'HandleVisibility','off');
                end
            catch
            end
        end

        xlim(ax,[0.6 max(occLevels)+0.4]);
        set(ax,'XTick',occLevels);
        grid(ax,'on');
        title(ax, sprintf('%s | WWR%s (n=%d)', cxv, wv, sum(idxPanel)), 'Interpreter','none');
        if r==numel(cxLevels)
            xlabel(ax,'Scene occurrence index');
        end
        if c==1
            ylabel(ax, char(metric));
        end

        if r==1 && c==numel(wwrLevels)
            legend(ax,'Location','best');
        end

        % remove toolbar
        try
            axtoolbar(ax, []);
        catch
        end
        try
            disableDefaultInteractivity(ax);
        catch
        end
    end
end

sgtitle(sprintf('Scene-occurrence trend (faceted) | %s | %s [%s]', string(analysisName), string(metric), string(tag)), 'Interpreter','none');

% annotate p-values if available
% NOTE: these p-values come from the global TrialIndex LMM (1..12), not the
% panel-local occurrence-index visualization.
try
    if nargin>=6 && ~isempty(lme)
        C = lme.Coefficients;
        rn = string(C.Name);
        ann = "";
        i1 = find(rn=="TrialIndex",1);
        if ~isempty(i1)
            ann = ann + sprintf('Global TrialIndex p=%.3g; ', C.pValue(i1));
        end
        i2 = find(contains(rn,"TrialIndex") & contains(rn,"Group"),1);
        if ~isempty(i2)
            ann = ann + sprintf('Group×TrialIndex p=%.3g', C.pValue(i2));
        end
        if strlength(ann)>0
            annotation(fig,'textbox',[0.01 0.01 0.98 0.05],'String',ann,'EdgeColor','none','HorizontalAlignment','center');
        end
    end
catch
end

fp = fullfile(fp_out, pipeline.sanitize_filename(sprintf('sceneocc_trend_facets_%s_%s_%s.png', tag, analysisName, metric)));
pipeline.export_figure_png(fig, fp, 300);
try; close(fig); catch; end
end
