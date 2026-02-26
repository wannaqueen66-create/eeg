function fp = plot_trialindex_trend(Tready, analysisName, metric, tag, fp_out, lme)
%PLOT_TRIALINDEX_TREND Paper-friendly plot: mean EEG vs TrialIndex by Group.
%
% Tready columns: Group (categorical), TrialIndex (double), EEG (double)

if nargin < 3; metric = "EEG"; end
if nargin < 4; tag = "raw"; end
if nargin < 5 || strlength(string(fp_out))==0
    fp_out = pwd;
end

if ~exist(fp_out,'dir'); mkdir(fp_out); end

set(0,'DefaultFigureVisible','off');

% ensure types
try
    grp = string(Tready.Group);
catch
    grp = string(categorical(Tready.Group));
end
x = double(Tready.TrialIndex);
y = double(Tready.EEG);

% only High/Low
mask = (grp=="High" | grp=="Low") & ~isnan(x) & ~isnan(y);
grp = grp(mask); x = x(mask); y = y(mask);

if numel(y) < 20
    fp = "";
    return;
end

% summary by trial index
trials = (1:12)';
colors = struct('Low',[0.2 0.5 0.9], 'High',[0.9 0.4 0.2]);

fig = figure('Position',[100 100 980 420], 'Color','w');
ax = axes(fig); %#ok<LAXES>
hold(ax,'on');

for gname = ["Low","High"]
    idx = (grp==gname);
    if sum(idx) < 10
        continue;
    end
    mu = nan(numel(trials),1);
    sem = nan(numel(trials),1);
    for i=1:numel(trials)
        xi = trials(i);
        yy = y(idx & x==xi);
        if isempty(yy)
            continue;
        end
        mu(i) = mean(yy,'omitnan');
        sem(i) = std(yy,'omitnan')/sqrt(numel(yy));
    end

    % plot mean with error bars
    c = colors.(char(gname));
    errorbar(ax, trials, mu, sem, 'o-', 'LineWidth', 1.8, 'Color', c, ...
        'MarkerFaceColor', c, 'DisplayName', char(gname));

    % simple visual linear fit (not inference)
    try
        ok = ~isnan(mu);
        if sum(ok) >= 3
            p = polyfit(trials(ok), mu(ok), 1);
            yhat = polyval(p, trials);
            plot(ax, trials, yhat, '--', 'LineWidth', 1.3, 'Color', c, 'HandleVisibility','off');
        end
    catch
    end
end

xlim(ax,[0.6 12.4]);
set(ax,'XTick',1:12);
grid(ax,'on');
xlabel(ax,'TrialIndex (1–12)');
ylabel(ax, char(metric));

ttl = sprintf('TrialIndex trend | %s | %s [%s]', string(analysisName), string(metric), string(tag));

% annotate key p-values from LME if provided
ann = "";
try
    if nargin>=6 && ~isempty(lme)
        C = lme.Coefficients;
        rn = string(C.Name);
        % main TrialIndex
        i1 = find(rn=="TrialIndex",1);
        if ~isempty(i1)
            ann = ann + sprintf('TrialIndex: p=%.3g; ', C.pValue(i1));
        end
        % interaction term name depends on coding
        i2 = find(contains(rn,"TrialIndex") & contains(rn,"Group"),1);
        if ~isempty(i2)
            ann = ann + sprintf('Group×TrialIndex: p=%.3g', C.pValue(i2));
        end
    end
catch
end

if strlength(ann)>0
    title(ax, ttl + " | " + ann, 'Interpreter','none');
else
    title(ax, ttl, 'Interpreter','none');
end

% Export
fp = fullfile(fp_out, pipeline.sanitize_filename(sprintf('trialindex_trend_%s_%s_%s.png', tag, analysisName, metric)));
pipeline.export_figure_png(fig, fp, 300);
try; close(fig); catch; end
end
