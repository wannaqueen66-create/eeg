function out = plot_eye_eeg_links(AllSceneEye, fp_sum, cfg, tag)
%PLOT_EYE_EEG_LINKS Minimal figures for EEG x eye exploratory inspection.

if nargin < 4 || isempty(tag)
    tag = 'raw';
end
if nargin < 3
    cfg = struct();
end
out = struct();
if isempty(AllSceneEye) || ~istable(AllSceneEye) || height(AllSceneEye)==0
    return;
end

fp_fig = fullfile(fp_sum, 'analysis-eye', char(tag), 'figures');
if ~exist(fp_fig,'dir'); mkdir(fp_fig); end

EEG = {'F_theta','O_alpha','O_beta'};
EYE = {'eye_view_blink_rate_per_min','eye_view_sacc_rate_per_min','eye_tracking_ratio','eye_mean_pupil_mm'};
try
    if isfield(cfg,'eye_analysis_eeg_metrics') && ~isempty(cfg.eye_analysis_eeg_metrics)
        EEG = cellstr(string(cfg.eye_analysis_eeg_metrics));
    end
catch
end
try
    if isfield(cfg,'eye_analysis_eye_metrics') && ~isempty(cfg.eye_analysis_eye_metrics)
        EYE = cellstr(string(cfg.eye_analysis_eye_metrics));
    end
catch
end
EEG = EEG(ismember(EEG, AllSceneEye.Properties.VariableNames));
EYE = EYE(ismember(EYE, AllSceneEye.Properties.VariableNames));
if isempty(EEG) || isempty(EYE)
    return;
end

blue = [0.20 0.50 0.90];
orange = [0.90 0.40 0.20];

% Scatter plots with cleaner journal-like styling
for i = 1:numel(EEG)
    for j = 1:numel(EYE)
        x = double(AllSceneEye.(EYE{j}));
        y = double(AllSceneEye.(EEG{i}));
        keep = ~isnan(x) & ~isnan(y);
        if sum(keep) < 8
            continue;
        end

        fig = figure('Visible','off','Color','w','Position',[80 80 720 520]);
        ax = axes(fig); hold(ax,'on');
        scatter(ax, x(keep), y(keep), 26, blue, 'filled', 'MarkerFaceAlpha',0.55, 'MarkerEdgeColor','none');

        p = polyfit(x(keep), y(keep), 1);
        xx = linspace(min(x(keep)), max(x(keep)), 100);
        yy = polyval(p, xx);
        plot(ax, xx, yy, '-', 'Color', orange, 'LineWidth',1.8);

        try
            [r,pp] = corr(x(keep), y(keep), 'Type','Spearman', 'Rows','complete');
            ttl = sprintf('%s vs %s | r=%.2f, p=%.3g, n=%d [%s]', EEG{i}, EYE{j}, r, pp, sum(keep), tag);
        catch
            ttl = sprintf('%s vs %s | n=%d [%s]', EEG{i}, EYE{j}, sum(keep), tag);
        end
        title(ax, ttl, 'Interpreter','none', 'FontWeight','normal');
        xlabel(ax, strrep(EYE{j}, '_', '\_'), 'Interpreter','tex');
        ylabel(ax, EEG{i}, 'Interpreter','none');
        grid(ax,'on'); box(ax,'off');
        ax.GridAlpha = 0.18;
        ax.MinorGridAlpha = 0.10;
        ax.FontName = 'Helvetica';
        ax.FontSize = 11;
        ax.LineWidth = 1;

        fp = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('scatter_%s_vs_%s.png', EEG{i}, EYE{j})));
        pipeline.export_figure_png(fig, fp, cfg);
        close(fig);
    end
end

% Box/scatter by eye_qc_needs_review flag if present
if ismember('eye_qc_needs_review', AllSceneEye.Properties.VariableNames)
    g = double(AllSceneEye.eye_qc_needs_review);
    keepg = ~isnan(g);
    for i = 1:numel(EEG)
        if ~ismember(EEG{i}, AllSceneEye.Properties.VariableNames)
            continue;
        end
        y = double(AllSceneEye.(EEG{i}));
        keep = keepg & ~isnan(y);
        if sum(keep) < 8 || numel(unique(g(keep))) < 2
            continue;
        end

        fig = figure('Visible','off','Color','w','Position',[70 70 900 460]);
        ax = axes(fig); hold(ax,'on');
        cols = [blue; orange];

        hOk = scatter(ax, nan, nan, 26, cols(1,:), 'filled', 'MarkerFaceAlpha',0.55, 'DisplayName','QC ok');
        hRv = scatter(ax, nan, nan, 26, cols(2,:), 'filled', 'MarkerFaceAlpha',0.55, 'DisplayName','Needs review');

        for gi = 0:1
            idx = keep & g==gi;
            yy = y(idx);
            if isempty(yy); continue; end
            x0 = gi + 1;
            try
                boxchart(ax, repmat(x0, numel(yy),1), yy, 'BoxWidth',0.28, 'MarkerStyle','none', ...
                    'BoxFaceColor', cols(gi+1,:), 'BoxFaceAlpha',0.22, 'HandleVisibility','off');
            catch
                q = quantile(yy,[0.25 0.5 0.75]);
                plot(ax,[x0-0.10 x0+0.10],[q(2) q(2)],'-','Color',cols(gi+1,:),'LineWidth',2, 'HandleVisibility','off');
            end
            xx = repmat(x0, numel(yy), 1) + (rand(size(yy))-0.5)*0.12;
            scatter(ax, xx, yy, 18, cols(gi+1,:), 'filled', 'MarkerFaceAlpha',0.50, 'HandleVisibility','off');
        end

        set(ax, 'XTick', [1 2], 'XTickLabel', {'QC ok','Needs review'});
        ylabel(ax, EEG{i}, 'Interpreter','none');
        xlabel(ax, 'Eye QC grouping');
        title(ax, sprintf('Eye-QC grouped box/scatter | %s [%s]', EEG{i}, tag), 'Interpreter','none', 'FontWeight','normal');
        grid(ax,'on'); box(ax,'off');
        ax.GridAlpha = 0.18;
        ax.FontName = 'Helvetica';
        ax.FontSize = 11;
        ax.LineWidth = 1;
        legend(ax, [hOk hRv], 'Location','best');

        fp = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('box_%s_by_eye_qc_review.png', EEG{i})));
        pipeline.export_figure_png(fig, fp, cfg);
        close(fig);
    end
end

out.fig_dir = fp_fig;
end
