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

fig = figure('Position',[100 100 780 420], 'Color','w');
hold on;
for i=1:n
    plot([1 2], [y1(i) y2(i)], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
end
scatter(ones(n,1)*1, y1, 30, 'filled', 'MarkerFaceColor',[0.2 0.5 0.9]);
scatter(ones(n,1)*2, y2, 30, 'filled', 'MarkerFaceColor',[0.9 0.4 0.2]);

mu1 = mean(y1,'omitnan'); mu2 = mean(y2,'omitnan');
se1 = std(y1,'omitnan')/sqrt(n); se2 = std(y2,'omitnan')/sqrt(n);
errorbar(1, mu1, se1, 'k', 'LineWidth', 2);
errorbar(2, mu2, se2, 'k', 'LineWidth', 2);

xlim([0.6 2.4]);
set(gca,'XTick',[1 2],'XTickLabel',{label1,label2});
grid on;
title(ttl, 'Interpreter','none');
ylabel(ylab);

pipeline.export_figure_png(fig, fp_out, 300);
try; close(fig); catch; end
end
