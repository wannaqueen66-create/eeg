function plot_paired_summary(y1, y2, label1, label2, ttl, ylab, fp_out)
%PLOT_PAIRED_SUMMARY Paper-friendly paired plot with subject lines + mean±SEM.

set(0,'DefaultFigureVisible','off');

y1 = double(y1(:));
y2 = double(y2(:));
mask = ~isnan(y1) & ~isnan(y2);
y1 = y1(mask);
y2 = y2(mask);

n = numel(y1);
if n < 3
    return;
end

fig = figure('Position',[100 100 820 430], 'Color','w');
ax = axes(fig); hold(ax,'on');
for i=1:n
    plot(ax, [1 2], [y1(i) y2(i)], '-', 'Color', [0.80 0.80 0.80], 'LineWidth', 0.9);
end
scatter(ax, ones(n,1)*1, y1, 30, 'filled', 'MarkerFaceColor',[0.2 0.5 0.9], 'MarkerFaceAlpha',0.65);
scatter(ax, ones(n,1)*2, y2, 30, 'filled', 'MarkerFaceColor',[0.9 0.4 0.2], 'MarkerFaceAlpha',0.65);

mu1 = mean(y1,'omitnan'); mu2 = mean(y2,'omitnan');
se1 = std(y1,'omitnan')/sqrt(n); se2 = std(y2,'omitnan')/sqrt(n);
errorbar(ax, 1, mu1, se1, 'o', 'Color',[0.2 0.5 0.9], 'MarkerFaceColor',[0.2 0.5 0.9], 'LineWidth', 1.5);
errorbar(ax, 2, mu2, se2, 'o', 'Color',[0.9 0.4 0.2], 'MarkerFaceColor',[0.9 0.4 0.2], 'LineWidth', 1.5);

xlim(ax,[0.6 2.4]);
set(ax,'XTick',[1 2],'XTickLabel',{label1,label2});
ax.FontName = 'Times New Roman';
ax.FontSize = 11;
ax.LineWidth = 0.8;
box(ax,'off');
grid(ax,'on');
ax.GridAlpha = 0.08;
ax.GridColor = [0 0 0];
title(ax, ttl, 'Interpreter','none', 'FontWeight','normal');
ylabel(ax, ylab);

pipeline.export_figure_png(fig, fp_out, 300);
try; close(fig); catch; end
end
