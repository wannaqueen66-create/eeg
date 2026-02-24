function plot_group_scene_by_factors(AllScene, fp_sum, cfg, tag)
%PLOT_GROUP_SCENE_BY_FACTORS Journal-friendly grouped plots by WWR and Complexity.
%
% Outputs (under summary/fig):
%   B1) Complexity facets + WWR ordering
%   B2) WWR × Complexity interaction plot
%
% For each factor (Experience/SportFreq) and each main metric.

if nargin < 4; tag = ""; end
if nargin < 3; cfg = struct(); end

if isempty(AllScene) || ~istable(AllScene) || height(AllScene)==0
    return;
end

needCols = {'subject_id','Complexity'};
for i=1:numel(needCols)
    if ~ismember(needCols{i}, AllScene.Properties.VariableNames)
        warning('plot_group_scene_by_factors: missing %s', needCols{i});
        return;
    end
end

% WWR column name variants
wcol = '';
if ismember('WWR', AllScene.Properties.VariableNames)
    wcol = 'WWR';
elseif ismember('WWRlevel', AllScene.Properties.VariableNames)
    wcol = 'WWRlevel';
end
if isempty(wcol)
    warning('plot_group_scene_by_factors: missing WWR');
    return;
end

fp_fig = fullfile(fp_sum, 'fig');
if ~exist(fp_fig,'dir'); mkdir(fp_fig); end

metrics = {"O_alpha","O_theta","O_beta","F_theta"};
try
    if isfield(cfg,'group_plot_metrics') && ~isempty(cfg.group_plot_metrics)
        metrics = cellstr(string(cfg.group_plot_metrics));
    elseif isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
catch
end

cx01 = local_complexity01(AllScene);
WWR = AllScene.(wcol);
if iscell(WWR); WWR = string(WWR); end
if isstring(WWR)
    WWRv = double(WWR);
else
    WWRv = double(WWR);
end

factors = {"Experience","SportFreq"};

for fi=1:numel(factors)
    fac = factors{fi};
    if ~ismember(fac, AllScene.Properties.VariableNames)
        continue;
    end

    grp = strtrim(string(AllScene.(fac)));
    grp(lower(grp)=="high") = "High";
    grp(lower(grp)=="low") = "Low";

    keep0 = (grp=="High" | grp=="Low") & ~isnan(cx01) & ~isnan(WWRv);
    if sum(keep0) < 6
        continue;
    end

    A = AllScene(keep0,:);
    A.group_attach = grp(keep0);
    A.cx01_attach = cx01(keep0);
    A.WWR_attach = WWRv(keep0);

    sid = string(A.subject_id);

    for mi=1:numel(metrics)
        met = string(metrics{mi});
        if ~ismember(met, A.Properties.VariableNames)
            continue;
        end

        y = A.(met);
        if ~isnumeric(y)
            try; y = double(y); catch; continue; end
        end

        % subject-level mean within each (group, cx, WWR)
        [G1, sid_u, grp_u, cx_u, w_u] = findgroups(sid, A.group_attach, A.cx01_attach, A.WWR_attach);
        mu = splitapply(@(x) mean(x,'omitnan'), y, G1);
        Tsub = table(sid_u, grp_u, cx_u, w_u, mu, 'VariableNames', {'subject_id','group','cx','WWR','mu'});

        % group mean/sem
        [G2, gg, cc, ww] = findgroups(Tsub.group, Tsub.cx, Tsub.WWR);
        gmu  = splitapply(@(x) mean(x,'omitnan'), Tsub.mu, G2);
        gsem = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), Tsub.mu, G2);
        gn   = splitapply(@(x) sum(~isnan(x)), Tsub.mu, G2);
        Tsum = table(gg, cc, ww, gmu, gsem, gn, 'VariableNames', {'group','cx','WWR','mean','sem','n'});

        wLevels = sort(unique(Tsum.WWR));
        if isempty(wLevels); continue; end

        %% B1: Complexity facets + WWR ordering
        fig = figure('Name', sprintf('%s %s B1', fac, met), 'Position', [100 100 1200 420]);
        t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        title(t, sprintf('%s: %s by WWR (faceted by Complexity)%s', fac, met, local_tag(tag)), 'Interpreter','none');

        for cxi=0:1
            nexttile;
            hold on; grid on;
            cols = lines(2);
            groups_order = ["Low","High"];
            for gi=1:2
                glab = groups_order(gi);
                Ti = Tsum(Tsum.group==glab & Tsum.cx==cxi,:);
                if isempty(Ti); continue; end
                [~,o] = sort(Ti.WWR); Ti = Ti(o,:);
                errorbar(Ti.WWR, Ti.mean, Ti.sem, 'o-', 'LineWidth', 1.8, 'Color', cols(gi,:), ...
                    'MarkerFaceColor', cols(gi,:), 'DisplayName', sprintf('%s (mean n≈%.1f)', glab, mean(Ti.n,'omitnan')));
            end
            set(gca,'XTick',wLevels);
            xlabel('WWR');
            ylabel(char(met), 'Interpreter','none');
            if cxi==0
                title('ComplexityLow');
            else
                title('ComplexityHigh');
            end
            legend('Location','best');
        end

        fn = sprintf('group_B1_WWRfacets_%s_%s', lower(fac), char(met));
        if strlength(string(tag))>0; fn = fn + "_" + string(tag); end
        saveas(fig, fullfile(fp_fig, pipeline.sanitize_filename(fn + ".png")));
        try; close(fig); catch; end

        %% B2: WWR × Complexity interaction plot
        fig = figure('Name', sprintf('%s %s B2', fac, met), 'Position', [100 100 900 520]);
        hold on; grid on;
        cols = lines(4);
        % line styles: group + complexity
        combos = {"Low",0; "Low",1; "High",0; "High",1};
        for ci=1:size(combos,1)
            glab = combos{ci,1};
            cxi  = combos{ci,2};
            Ti = Tsum(Tsum.group==glab & Tsum.cx==cxi,:);
            if isempty(Ti); continue; end
            [~,o] = sort(Ti.WWR); Ti = Ti(o,:);
            if cxi==0
                ls = '-';
                cxlab = 'CxLow';
            else
                ls = '--';
                cxlab = 'CxHigh';
            end
            errorbar(Ti.WWR, Ti.mean, Ti.sem, 'o', 'LineWidth', 1.5, 'Color', cols(ci,:), 'MarkerFaceColor', cols(ci,:));
            plot(Ti.WWR, Ti.mean, ls, 'LineWidth', 1.8, 'Color', cols(ci,:), 'DisplayName', sprintf('%s-%s', glab, cxlab));
        end
        set(gca,'XTick',wLevels);
        xlabel('WWR');
        ylabel(char(met), 'Interpreter','none');
        title(sprintf('%s: %s WWR×Complexity interaction%s', fac, met, local_tag(tag)), 'Interpreter','none');
        legend('Location','best');

        fn = sprintf('group_B2_WWRxComplexity_%s_%s', lower(fac), char(met));
        if strlength(string(tag))>0; fn = fn + "_" + string(tag); end
        saveas(fig, fullfile(fp_fig, pipeline.sanitize_filename(fn + ".png")));
        try; close(fig); catch; end

    end
end

end

function cx01 = local_complexity01(T)
% normalize Complexity to 0/1
cx01 = nan(height(T),1);
if ~ismember('Complexity', T.Properties.VariableNames)
    return;
end
cx_raw = T.Complexity;
if iscell(cx_raw); cx_raw = string(cx_raw); end
if isstring(cx_raw)
    cx = lower(strtrim(cx_raw));
    cx01(cx=="low" | cx=="0") = 0;
    cx01(cx=="high" | cx=="1") = 1;
else
    cx01 = double(cx_raw);
end
end

function s = local_tag(tag)
if strlength(string(tag))>0
    s = " [" + string(tag) + "]";
else
    s = "";
end
end
