function plot_paper_figures(fp_sum, cfg)
%PLOT_PAPER_FIGURES Create journal-ready multi-panel figures under summary/paper_fig/.
%
% Uses summary merged CSVs if present:
%   - all_subjects_scene_level.csv
%   - all_subjects_pairs_check.csv
%
% Output:
%   summary/paper_fig/

if nargin < 2
    cfg = struct();
end

fp_scene = fullfile(fp_sum, 'all_subjects_scene_level.csv');
fp_pairs = fullfile(fp_sum, 'all_subjects_pairs_check.csv');

if ~exist(fp_scene,'file')
    warning('plot_paper_figures: missing %s', fp_scene);
    return;
end

S = readtable(fp_scene);

P = table();
if exist(fp_pairs,'file')
    try
        P = readtable(fp_pairs);
    catch
    end
end

fp_fig = fullfile(fp_sum, 'paper_fig');
if ~exist(fp_fig,'dir'); mkdir(fp_fig); end

% Style
try
    set(0,'DefaultAxesFontName','Arial');
    set(0,'DefaultTextFontName','Arial');
    set(0,'DefaultAxesFontSize',11);
end

paper_dpi = 300;
try
    if isfield(cfg,'paper_dpi') && ~isempty(cfg.paper_dpi)
        paper_dpi = double(cfg.paper_dpi);
    end
catch
end

metrics = {"O_alpha","O_theta","O_beta","F_theta"};
try
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
catch
end

% Helper: normalize Complexity to 0/1
cx01 = local_complexity01(S);
S.cx01 = cx01;

% Ensure groups exist
factors = ["Experience","SportFreq"];
for fi=1:numel(factors)
    fac = factors(fi);
    if ~ismember(fac, S.Properties.VariableNames)
        continue;
    end

    % Build subject means per complexity
    sid = string(S.subject_id);
    grp = strtrim(string(S.(fac)));
    grp(lower(grp)=="high") = "High";
    grp(lower(grp)=="low") = "Low";

    keep = ~isnan(S.cx01) & (grp=="High" | grp=="Low");
    if sum(keep) < 10
        continue;
    end

    % Multi-panel: mean by complexity
    fig = figure('Name', sprintf('Paper %s by Complexity', fac), 'Position',[100 100 1100 760]);
    t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(t, sprintf('%s groups × Complexity (subject means)', fac), 'Interpreter','none');

    for mi=1:numel(metrics)
        met = char(metrics{mi});
        if ~ismember(met, S.Properties.VariableNames)
            continue;
        end

        nexttile;
        y = S.(met);
        if ~isnumeric(y)
            try; y = double(y); catch; continue; end
        end

        % subject-level mean for each complexity
        [G1, sid_u, grp_u, cx_u] = findgroups(sid(keep), grp(keep), S.cx01(keep));
        mu = splitapply(@(x) mean(x,'omitnan'), y(keep), G1);
        Tsub = table(sid_u, grp_u, cx_u, mu, 'VariableNames', {'subject_id','group','cx','mu'});

        local_plot_by_complexity(Tsub, met);
    end

    local_export(fig, fullfile(fp_fig, sprintf('Fig_%s_by_Complexity', lower(fac))), paper_dpi);

    % Multi-panel: delta complexity (High-Low) by group
    fig = figure('Name', sprintf('Paper %s Delta Complexity', fac), 'Position',[100 100 1100 760]);
    t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(t, sprintf('%s groups: Complexity effect (High-Low) per subject', fac), 'Interpreter','none');

    for mi=1:numel(metrics)
        met = char(metrics{mi});
        if ~ismember(met, S.Properties.VariableNames)
            continue;
        end

        nexttile;
        y = S.(met);
        if ~isnumeric(y)
            try; y = double(y); catch; continue; end
        end

        [G1, sid_u, grp_u, cx_u] = findgroups(sid(keep), grp(keep), S.cx01(keep));
        mu = splitapply(@(x) mean(x,'omitnan'), y(keep), G1);
        Tsub = table(sid_u, grp_u, cx_u, mu, 'VariableNames', {'subject_id','group','cx','mu'});

        local_plot_delta_complexity(Tsub, met);
    end

    local_export(fig, fullfile(fp_fig, sprintf('Fig_%s_delta_Complexity', lower(fac))), paper_dpi);

    % Recovery paper plots (delta_O_alpha) if available
    if ~isempty(P)
        try
            cxp = local_complexity01(P);
            if ismember('delta_O_alpha', P.Properties.VariableNames)
                fig = figure('Name', sprintf('Paper %s Recovery', fac), 'Position',[100 100 1100 420]);
                tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

                % By complexity
                nexttile;
                local_plot_recovery_by_complexity(P, fac, cxp);

                % Delta complexity
                nexttile;
                local_plot_recovery_delta_complexity(P, fac, cxp);

                local_export(fig, fullfile(fp_fig, sprintf('Fig_%s_recovery_delta_O_alpha', lower(fac))), paper_dpi);
            end
        catch
        end
    end
end

end

function cx01 = local_complexity01(T)
if ~ismember('Complexity', T.Properties.VariableNames)
    cx01 = nan(height(T),1);
    return;
end
cx_raw = T.Complexity;
if iscell(cx_raw); cx_raw = string(cx_raw); end
if isstring(cx_raw)
    cx = lower(strtrim(cx_raw));
    cx01 = nan(height(T),1);
    cx01(cx=="low" | cx=="0") = 0;
    cx01(cx=="high" | cx=="1") = 1;
else
    cx01 = double(cx_raw);
end
end

function local_plot_by_complexity(Tsub, met)
% Tsub: subject_id, group(High/Low), cx(0/1), mu
hold on; grid on;
groups = ["Low","High"]; cols = lines(2);
for gi=1:2
    g = groups(gi);
    A = Tsub(Tsub.group==g,:);
    if isempty(A), continue; end
    [G2, cx_u] = findgroups(A.cx);
    m = splitapply(@(x) mean(x,'omitnan'), A.mu, G2);
    s = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), A.mu, G2);
    [cx_u, o] = sort(cx_u);
    m=m(o); s=s(o);
    errorbar(cx_u, m, s, 'o-', 'LineWidth', 1.8, 'Color', cols(gi,:), 'MarkerFaceColor', cols(gi,:), 'DisplayName', g);
end
set(gca,'XTick',[0 1],'XTickLabel',{'ComplexityLow','ComplexityHigh'});
xlabel('Complexity');
ylabel(met, 'Interpreter','none');
title(met, 'Interpreter','none');
legend('Location','best');
end

function local_plot_delta_complexity(Tsub, met)
% Delta per subject
% Compute subject delta = High - Low
[Gs, sid_u, grp_u] = findgroups(Tsub.subject_id, Tsub.group);
del = splitapply(@(cx, mu) local_delta(cx, mu), Tsub.cx, Tsub.mu, Gs);
D = table(sid_u, grp_u, del, 'VariableNames', {'subject_id','group','delta'});

hold on; grid on;
% group mean/sem
[G2, g_u] = findgroups(D.group);
m = splitapply(@(x) mean(x,'omitnan'), D.delta, G2);
s = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), D.delta, G2);
order = ["Low";"High"];
x = 1:2;
mm = nan(2,1); ss = nan(2,1);
for i=1:2
    k = find(g_u==order(i));
    if ~isempty(k)
        mm(i)=m(k);
        ss(i)=s(k);
    end
end
bar(x, mm, 'FaceColor',[0.8 0.85 0.95]);
errorbar(x, mm, ss, 'k.', 'LineWidth', 1.5);
set(gca,'XTick',x,'XTickLabel',order);
ylabel(sprintf('%s (ComplexityHigh - ComplexityLow)', met), 'Interpreter','none');
title(met, 'Interpreter','none');
end

function d = local_delta(cx, mu)
lo = mu(cx==0);
hi = mu(cx==1);
if isempty(lo) || isempty(hi)
    d = NaN;
else
    d = mean(hi,'omitnan') - mean(lo,'omitnan');
end
end

function local_plot_recovery_by_complexity(P, fac, cx01)
% P includes subject_id, delta_O_alpha, fac, Complexity
sid = string(P.subject_id);
grp = strtrim(string(P.(fac)));
grp(lower(grp)=="high") = "High";
grp(lower(grp)=="low") = "Low";
keep = ~isnan(cx01) & (grp=="High" | grp=="Low") & ~isnan(P.delta_O_alpha);

[Gs, sid_u, grp_u, cx_u] = findgroups(sid(keep), grp(keep), cx01(keep));
mu = splitapply(@(x) mean(x,'omitnan'), P.delta_O_alpha(keep), Gs);
Tsub = table(sid_u, grp_u, cx_u, mu, 'VariableNames', {'subject_id','group','cx','mu'});

hold on; grid on;
local_plot_by_complexity(Tsub, 'delta_O_alpha');
title(sprintf('%s recovery by Complexity', fac), 'Interpreter','none');
end

function local_plot_recovery_delta_complexity(P, fac, cx01)
% Delta complexity of recovery
sid = string(P.subject_id);
grp = strtrim(string(P.(fac)));
grp(lower(grp)=="high") = "High";
grp(lower(grp)=="low") = "Low";
keep = ~isnan(cx01) & (grp=="High" | grp=="Low") & ~isnan(P.delta_O_alpha);

[Gs, sid_u, grp_u, cx_u] = findgroups(sid(keep), grp(keep), cx01(keep));
mu = splitapply(@(x) mean(x,'omitnan'), P.delta_O_alpha(keep), Gs);
Tsub = table(sid_u, grp_u, cx_u, mu, 'VariableNames', {'subject_id','group','cx','mu'});

hold on; grid on;
local_plot_delta_complexity(Tsub, 'delta_O_alpha');
title(sprintf('%s recovery ΔComplexity', fac), 'Interpreter','none');
end

function local_export(fig, fp_base, dpi)
% Save both PNG and PDF
try
    set(fig,'Color','w');
catch
end
try
    print(fig, fp_base + ".png", '-dpng', ['-r' num2str(dpi)]);
catch
    try
        saveas(fig, fp_base + ".png");
    catch
    end
end
try
    print(fig, fp_base + ".pdf", '-dpdf');
catch
end
try; close(fig); catch; end
end
