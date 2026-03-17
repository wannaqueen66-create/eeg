function plot_group_summaries(AllScene, fp_batch, cfg, tag)
%PLOT_GROUP_SUMMARIES Generate group-level PNG summaries under summary figures.
%
% Uses merged scene-level table (typically summary/all_subjects_scene_level.csv)
% and produces simple group comparison plots:
%   - Experience (High vs Low) by Complexity (Low vs High)
%   - SportFreq (High vs Low) by Complexity (Low vs High)
%
% Plots are based on subject-level means (within-subject averaging across scenes
% of the same complexity), then group means/SEM across subjects.

if nargin < 3; cfg = struct(); end
if nargin < 4; tag = 'raw'; end

if isempty(AllScene) || ~istable(AllScene) || height(AllScene)==0
    return;
end

if ~ismember('subject_id', AllScene.Properties.VariableNames)
    warning('plot_group_summaries: missing subject_id');
    return;
end

% output dir
fp_fig = pipeline.get_fig_dir(fp_batch, cfg, 'group', tag);

% metrics list (configurable)
% Journal-oriented defaults (scene effects in occipital bands + a frontal control)
default_metrics = {"O_alpha","O_theta","O_beta","F_theta"};
metrics = default_metrics;
try
    if isfield(cfg,'group_plot_metrics') && ~isempty(cfg.group_plot_metrics)
        m = cfg.group_plot_metrics;
        if ischar(m) || isstring(m)
            metrics = cellstr(string(m));
        elseif iscell(m)
            metrics = cellstr(string(m));
        end
    end
catch
end

% Complexity column required
if ~ismember('Complexity', AllScene.Properties.VariableNames)
    warning('plot_group_summaries: Complexity column missing. Skipping group plots.');
    return;
end

% Normalize Complexity to Low/High label and numeric 0/1
cx_raw = AllScene.Complexity;
if iscell(cx_raw)
    cx_raw = string(cx_raw);
end
if isstring(cx_raw)
    cx = lower(strtrim(cx_raw));
    cx_low = cx=="low" | cx=="0";
    cx_high = cx=="high" | cx=="1";
    cx01 = nan(height(AllScene),1);
    cx01(cx_low) = 0;
    cx01(cx_high) = 1;
else
    cx01 = double(cx_raw);
end

% Factors (between-subject) to plot
factors = {"Experience","SportFreq"};

for fi = 1:numel(factors)
    fac = factors{fi};
    if ~ismember(fac, AllScene.Properties.VariableNames)
        continue;
    end

    group_raw = AllScene.(fac);
    if iscell(group_raw)
        group_raw = string(group_raw);
    end
    group = string(group_raw);
    group = strtrim(group);

    % keep only High/Low rows
    g_low = lower(group)=="low";
    g_high = lower(group)=="high";

    keep = ~isnan(cx01) & (g_low | g_high);
    if sum(keep) < 6
        continue;
    end

    A = AllScene(keep, :);
    % NOTE: table variable names must start with a letter.
    A.cx01_attach = cx01(keep);
    A.group_attach = group(keep);

    % Standardize group casing to High/Low
    A.group_attach(lower(A.group_attach)=="high") = "High";
    A.group_attach(lower(A.group_attach)=="low") = "Low";

    for mi = 1:numel(metrics)
        met = char(metrics{mi});
        if ~ismember(met, A.Properties.VariableNames)
            continue;
        end

        y = A.(met);
        if ~isnumeric(y)
            try; y = double(y); catch; continue; end
        end

        % subject-level means by group x complexity
        sid = string(A.subject_id);
        g = string(A.group_attach);
        cx = double(A.cx01_attach);

        [Gid, sid_u, g_u, cx_u] = findgroups(sid, g, cx);
        mu = splitapply(@(x) mean(x,'omitnan'), y, Gid);
        S = table(sid_u, g_u, cx_u, mu, 'VariableNames', {'subject_id','group','cx','mu'});

        % ensure both complexity levels exist
        if ~any(S.cx==0) || ~any(S.cx==1)
            continue;
        end

        % group-level mean/sem across subjects
        [G2, gg, cc] = findgroups(S.group, S.cx);
        gmu = splitapply(@(x) mean(x,'omitnan'), S.mu, G2);
        gsem = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), S.mu, G2);
        gn = splitapply(@(x) sum(~isnan(x)), S.mu, G2);

        Tsum = table(gg, cc, gmu, gsem, gn, 'VariableNames', {'group','cx','mean','sem','n'});

        % ---- Plot 1: mean by complexity ----
        fig = figure('Name', sprintf('%s %s by Complexity', fac, met), 'Position', [100 100 900 420]);
        hold on;
        grid on;

        % sort groups Low then High
        groups_order = ["Low","High"];
        colors = lines(2);

        for gi = 1:numel(groups_order)
            glab = groups_order(gi);
            Ti = Tsum(Tsum.group==glab, :);
            if isempty(Ti)
                continue;
            end
            % order by cx
            [~,o] = sort(Ti.cx);
            Ti = Ti(o,:);
            x = Ti.cx;
            errorbar(x, Ti.mean, Ti.sem, 'o-', 'LineWidth', 1.8, 'Color', colors(gi,:), ...
                'MarkerFaceColor', colors(gi,:), 'DisplayName', sprintf('%s (n=%d/%d)', glab, sum(Ti.n), numel(unique(S.subject_id(S.group==glab)))) );
        end

        set(gca,'XTick',[0 1],'XTickLabel',{'ComplexityLow','ComplexityHigh'});
        xlabel('Complexity');
        ylabel(met);
        title(sprintf('%s: %s by Complexity (subject means)', fac, met), 'Interpreter','none');
        legend('Location','northeastoutside');

        f1 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('group_%s_%s_by_Complexity.png', lower(fac), met)));
        saveas(fig, f1);
        try; close(fig); catch; end

        % ---- Plot 2: delta (ComplexityHigh - ComplexityLow) by group ----
        % compute delta per subject
        [G3, sid3, g3] = findgroups(S.subject_id, S.group);
        % For each subject, pick cx==0 and cx==1 means
        del = splitapply(@(cxv, muv) local_delta(cxv, muv), S.cx, S.mu, G3);
        D = table(sid3, g3, del, 'VariableNames', {'subject_id','group','delta'});

        % group stats
        [G4, g4] = findgroups(D.group);
        dmu = splitapply(@(x) mean(x,'omitnan'), D.delta, G4);
        dsem = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), D.delta, G4);
        dn = splitapply(@(x) sum(~isnan(x)), D.delta, G4);
        dT = table(g4, dmu, dsem, dn, 'VariableNames', {'group','mean','sem','n'});

        % enforce Low/High order
        dT = sortrows(dT, 'group');
        if any(dT.group=="Low") && any(dT.group=="High")
            dT = [dT(dT.group=="Low",:); dT(dT.group=="High",:)];
        end

        fig = figure('Name', sprintf('%s %s DeltaCx', fac, met), 'Position', [100 100 700 420]);
        hold on; grid on;
        x = 1:height(dT);
        bar(x, dT.mean, 'FaceColor',[0.7 0.7 0.9]);
        errorbar(x, dT.mean, dT.sem, 'k.', 'LineWidth', 1.5);
        set(gca, 'XTick', x, 'XTickLabel', dT.group);
        ylabel(sprintf('%s (ComplexityHigh - ComplexityLow)', met), 'Interpreter','none');
        title(sprintf('%s: %s Complexity effect by group', fac, met), 'Interpreter','none');

        for ii=1:height(dT)
            text(ii, dT.mean(ii) + max(0.02, abs(dT.mean(ii))*0.04 + dT.sem(ii)), sprintf('n=%d', dT.n(ii)), 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
        end

        f2 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('group_%s_%s_delta_Complexity.png', lower(fac), met)));
        saveas(fig, f2);
        try; close(fig); catch; end

    end
end

end

function d = local_delta(cx, mu)
% compute High(1) - Low(0) within a subject, robust to missing
cx = double(cx);
mu = double(mu);
lo = mu(cx==0);
hi = mu(cx==1);
if isempty(lo) || isempty(hi)
    d = NaN;
else
    d = mean(hi,'omitnan') - mean(lo,'omitnan');
end
end
