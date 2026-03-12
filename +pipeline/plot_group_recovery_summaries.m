function plot_group_recovery_summaries(AllPairs, fp_batch, cfg, tag)
%PLOT_GROUP_RECOVERY_SUMMARIES Group-level plots for view→gray recovery deltas.
%
% Reads merged pairs_check table and generates PNGs under summary figures.
%
% Current focus: delta_O_alpha (and optionally delta_O_theta, delta_O_beta if present)
% stratified by:
%   - Experience (High/Low)
%   - SportFreq (High/Low)
% and optionally by Complexity (Low/High) if available.

if nargin < 3; cfg = struct(); end
if nargin < 4; tag = 'raw'; end

if isempty(AllPairs) || ~istable(AllPairs) || height(AllPairs)==0
    return;
end

if ~ismember('subject_id', AllPairs.Properties.VariableNames)
    warning('plot_group_recovery_summaries: missing subject_id');
    return;
end

fp_fig = pipeline.get_fig_dir(fp_batch, cfg, 'recovery', tag);

% metrics to plot (pair-level deltas)
default_metrics = {"delta_O_alpha"};
metrics = default_metrics;
try
    if isfield(cfg,'group_recovery_metrics') && ~isempty(cfg.group_recovery_metrics)
        m = cfg.group_recovery_metrics;
        if ischar(m) || isstring(m)
            metrics = cellstr(string(m));
        elseif iscell(m)
            metrics = cellstr(string(m));
        end
    else
        % if other deltas exist, include them
        extra = {"delta_O_theta","delta_O_beta"};
        for i=1:numel(extra)
            if ismember(extra{i}, AllPairs.Properties.VariableNames)
                metrics{end+1} = extra{i}; %#ok<AGROW>
            end
        end
    end
catch
end

% Complexity optional
hasCx = ismember('Complexity', AllPairs.Properties.VariableNames);
if hasCx
    cx_raw = AllPairs.Complexity;
    if iscell(cx_raw); cx_raw = string(cx_raw); end
    if isstring(cx_raw)
        cx = lower(strtrim(cx_raw));
        cx_low = cx=="low" | cx=="0";
        cx_high = cx=="high" | cx=="1";
        cx01 = nan(height(AllPairs),1);
        cx01(cx_low)=0; cx01(cx_high)=1;
    else
        cx01 = double(cx_raw);
    end
else
    cx01 = nan(height(AllPairs),1);
end

factors = {"Experience","SportFreq"};
for fi=1:numel(factors)
    fac = factors{fi};
    if ~ismember(fac, AllPairs.Properties.VariableNames)
        continue;
    end

    graw = AllPairs.(fac);
    if iscell(graw); graw = string(graw); end
    group = string(graw); group = strtrim(group);
    g_low = lower(group)=="low";
    g_high = lower(group)=="high";

    for mi=1:numel(metrics)
        met = char(metrics{mi});
        if ~ismember(met, AllPairs.Properties.VariableNames)
            continue;
        end

        y = AllPairs.(met);
        if ~isnumeric(y)
            try; y=double(y); catch; continue; end
        end

        keep = (g_low|g_high) & ~isnan(y);
        if hasCx
            keep = keep & ~isnan(cx01);
        end
        if sum(keep) < 6
            continue;
        end

        P = AllPairs(keep,:);
        sid = string(P.subject_id);
        g = string(group(keep));
        g(lower(g)=="high")="High";
        g(lower(g)=="low")="Low";

        if hasCx
            cx = cx01(keep);
            % subject-level mean delta by cx within subject
            [G1, sid_u, g_u, cx_u] = findgroups(sid, g, cx);
            mu = splitapply(@(x) mean(x,'omitnan'), y(keep), G1);
            S = table(sid_u, g_u, cx_u, mu, 'VariableNames', {'subject_id','group','cx','mu'});

            % group-level mean/sem
            [G2, gg, cc] = findgroups(S.group, S.cx);
            gmu = splitapply(@(x) mean(x,'omitnan'), S.mu, G2);
            gsem = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), S.mu, G2);
            gn = splitapply(@(x) sum(~isnan(x)), S.mu, G2);
            Tsum = table(gg, cc, gmu, gsem, gn, 'VariableNames', {'group','cx','mean','sem','n'});

            % plot by complexity
            fig = figure('Name', sprintf('%s %s by Complexity', fac, met), 'Position', [100 100 900 420]);
            hold on; grid on;
            groups_order = ["Low","High"];
            colors = lines(2);
            for gi=1:2
                glab = groups_order(gi);
                Ti = Tsum(Tsum.group==glab,:);
                if isempty(Ti), continue; end
                [~,o]=sort(Ti.cx); Ti=Ti(o,:);
                errorbar(Ti.cx, Ti.mean, Ti.sem, 'o-', 'LineWidth', 1.8, 'Color', colors(gi,:), ...
                    'MarkerFaceColor', colors(gi,:), 'DisplayName', glab);
            end
            set(gca,'XTick',[0 1],'XTickLabel',{'ComplexityLow','ComplexityHigh'});
            xlabel('Complexity');
            ylabel(met, 'Interpreter','none');
            title(sprintf('Recovery (%s): %s by Complexity (subject means)', fac, met), 'Interpreter','none');
            legend('Location','best');

            f1 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('group_%s_%s_by_Complexity.png', lower(fac), met)));
            saveas(fig, f1);
            try; close(fig); catch; end

            % delta of delta: ComplexityHigh - ComplexityLow per subject
            [G3, sid3, g3] = findgroups(S.subject_id, S.group);
            del = splitapply(@(cxv, muv) local_delta(cxv, muv), S.cx, S.mu, G3);
            D = table(sid3, g3, del, 'VariableNames', {'subject_id','group','delta'});

            [G4, g4] = findgroups(D.group);
            dmu = splitapply(@(x) mean(x,'omitnan'), D.delta, G4);
            dsem = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), D.delta, G4);
            dn = splitapply(@(x) sum(~isnan(x)), D.delta, G4);
            dT = table(g4, dmu, dsem, dn, 'VariableNames', {'group','mean','sem','n'});
            if any(dT.group=="Low") && any(dT.group=="High")
                dT = [dT(dT.group=="Low",:); dT(dT.group=="High",:)];
            end

            fig = figure('Name', sprintf('%s %s DeltaComplexity', fac, met), 'Position', [100 100 700 420]);
            hold on; grid on;
            x = 1:height(dT);
            bar(x, dT.mean, 'FaceColor',[0.85 0.75 0.75]);
            errorbar(x, dT.mean, dT.sem, 'k.', 'LineWidth', 1.5);
            set(gca,'XTick',x,'XTickLabel',dT.group);
            ylabel(sprintf('%s: (ComplexityHigh - ComplexityLow)', met), 'Interpreter','none');
            title(sprintf('Recovery complexity effect by %s', fac), 'Interpreter','none');
            for ii=1:height(dT)
                text(ii, dT.mean(ii), sprintf(' n=%d', dT.n(ii)), 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
            end
            f2 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('group_%s_%s_delta_Complexity.png', lower(fac), met)));
            saveas(fig, f2);
            try; close(fig); catch; end

        else
            % no complexity: subject-level mean delta by group
            [G1, sid_u, g_u] = findgroups(sid, g);
            mu = splitapply(@(x) mean(x,'omitnan'), y(keep), G1);
            S = table(sid_u, g_u, mu, 'VariableNames', {'subject_id','group','mu'});

            [G2, gg] = findgroups(S.group);
            gmu = splitapply(@(x) mean(x,'omitnan'), S.mu, G2);
            gsem = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), S.mu, G2);
            gn = splitapply(@(x) sum(~isnan(x)), S.mu, G2);
            Tsum = table(gg, gmu, gsem, gn, 'VariableNames', {'group','mean','sem','n'});

            fig = figure('Name', sprintf('%s %s', fac, met), 'Position', [100 100 700 420]);
            hold on; grid on;
            if any(Tsum.group=="Low") && any(Tsum.group=="High")
                Tsum = [Tsum(Tsum.group=="Low",:); Tsum(Tsum.group=="High",:)];
            end
            x = 1:height(Tsum);
            bar(x, Tsum.mean, 'FaceColor',[0.75 0.85 0.75]);
            errorbar(x, Tsum.mean, Tsum.sem, 'k.', 'LineWidth', 1.5);
            set(gca,'XTick',x,'XTickLabel',Tsum.group);
            ylabel(met,'Interpreter','none');
            title(sprintf('Recovery (%s): %s', fac, met), 'Interpreter','none');
            for ii=1:height(Tsum)
                text(ii, Tsum.mean(ii), sprintf(' n=%d', Tsum.n(ii)), 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
            end
            f1 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('group_%s_%s.png', lower(fac), met)));
            saveas(fig, f1);
            try; close(fig); catch; end
        end
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
