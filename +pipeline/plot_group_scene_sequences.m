function plot_group_scene_sequences(AllScene, fp_sum, cfg, tag)
%PLOT_GROUP_SCENE_SEQUENCES Group-level scene sequence plots (round1..12) under summary/fig.
%
% Plots main metrics across scene_id with x-axis labels using scene_name/SceneID.
% For each factor (Experience/SportFreq): plot group mean±SEM per scene.
%
% Inputs:
%   AllScene: merged scene-level table (raw or QC-filtered)
%   fp_sum : summary output directory
%   cfg    : config struct (uses group_plot_metrics or paper_metrics)
%   tag    : optional string appended to filename (e.g., 'qc')

if nargin < 4; tag = ""; end
if nargin < 3; cfg = struct(); end

if isempty(AllScene) || ~istable(AllScene) || height(AllScene)==0
    return;
end

if ~ismember('subject_id', AllScene.Properties.VariableNames)
    warning('plot_group_scene_sequences: missing subject_id');
    return;
end
if ~ismember('scene_id', AllScene.Properties.VariableNames)
    warning('plot_group_scene_sequences: missing scene_id');
    return;
end

fp_fig = fullfile(fp_sum, 'fig');
if ~exist(fp_fig, 'dir'); mkdir(fp_fig); end

% metrics list
metrics = {"O_alpha","O_theta","O_beta","F_theta"};
try
    if isfield(cfg,'group_plot_metrics') && ~isempty(cfg.group_plot_metrics)
        metrics = cellstr(string(cfg.group_plot_metrics));
    elseif isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
catch
end

% scene label source
sceneLabels = local_scene_labels(AllScene);

% Complexity shading (optional)
hasCx = ismember('Complexity', AllScene.Properties.VariableNames);
if hasCx
    cx01 = local_complexity01(AllScene);
else
    cx01 = nan(height(AllScene),1);
end

factors = {"Experience","SportFreq"};

for fi=1:numel(factors)
    fac = factors{fi};
    if ~ismember(fac, AllScene.Properties.VariableNames)
        continue;
    end

    group_raw = AllScene.(fac);
    if iscell(group_raw); group_raw = string(group_raw); end
    group = strtrim(string(group_raw));
    group(lower(group)=="high") = "High";
    group(lower(group)=="low") = "Low";

    keep0 = (group=="High" | group=="Low") & ~isnan(AllScene.scene_id);
    if sum(keep0) < 6
        continue;
    end

    A = AllScene(keep0,:);
    A.group_attach = group(keep0);
    A.scene_label_attach = sceneLabels(keep0);
    if hasCx
        A.cx01_attach = cx01(keep0);
    else
        A.cx01_attach = nan(height(A),1);
    end

    % unique scenes (sorted)
    scenes = unique(double(A.scene_id));
    scenes = sort(scenes(:));

    % per scene complexity label for shading (majority vote)
    sceneCx = nan(numel(scenes),1);
    if hasCx
        for si=1:numel(scenes)
            v = A.cx01_attach(double(A.scene_id)==scenes(si));
            sceneCx(si) = mode(v(~isnan(v)));
        end
    end

    % tick labels (scene_name/SceneID) + journal-friendly shortening
    tickLabels = strings(numel(scenes),1);
    for si=1:numel(scenes)
        idx = double(A.scene_id)==scenes(si);
        lab = "";
        try
            lab = string(A.scene_label_attach(find(idx,1,'first')));
        catch
        end
        if strlength(lab)==0
            lab = "scene" + sprintf('%02d', scenes(si));
        end
        tickLabels(si) = local_shorten_label(lab);
    end

    for mi=1:numel(metrics)
        met = string(metrics{mi});
        if ~ismember(met, A.Properties.VariableNames)
            continue;
        end

        y = A.(met);
        if ~isnumeric(y)
            try; y = double(y); catch; continue; end
        end

        % subject-level mean per scene (in case duplicates)
        sid = string(A.subject_id);
        grp = string(A.group_attach);
        sc  = double(A.scene_id);

        [G1, sid_u, grp_u, sc_u] = findgroups(sid, grp, sc);
        mu = splitapply(@(x) mean(x,'omitnan'), y, G1);
        Tsub = table(sid_u, grp_u, sc_u, mu, 'VariableNames', {'subject_id','group','scene_id','mu'});

        % group mean/sem per scene
        [G2, gg, ss] = findgroups(Tsub.group, Tsub.scene_id);
        gmu  = splitapply(@(x) mean(x,'omitnan'), Tsub.mu, G2);
        gsem = splitapply(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))), Tsub.mu, G2);
        gn   = splitapply(@(x) sum(~isnan(x)), Tsub.mu, G2);
        Tsum = table(gg, ss, gmu, gsem, gn, 'VariableNames', {'group','scene_id','mean','sem','n'});

        fig = figure('Name', sprintf('%s %s scene sequence', fac, met), 'Position', [100 100 1200 420]);
        hold on; grid on;

        % optional background shading by Complexity
        try
            if hasCx && all(ismember([0 1], unique(sceneCx(~isnan(sceneCx)))'))
                yl = [min(Tsum.mean - Tsum.sem, [], 'omitnan'), max(Tsum.mean + Tsum.sem, [], 'omitnan')];
                if any(isnan(yl)); yl = ylim; end
                pad = 0.08 * (yl(2)-yl(1) + eps);
                yl = [yl(1)-pad, yl(2)+pad];
                ylim(yl);
                for si=1:numel(scenes)
                    if isnan(sceneCx(si)); continue; end
                    if sceneCx(si)==1
                        c = [1.0 0.95 0.95]; % high complexity
                    else
                        c = [0.95 0.97 1.0]; % low complexity
                    end
                    x0 = si-0.5; x1 = si+0.5;
                    patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], c, 'EdgeColor','none', 'FaceAlpha', 0.35);
                end
            end
        catch
        end

        groups_order = ["Low","High"];
        colors = lines(2);
        for gi=1:2
            glab = groups_order(gi);
            Ti = Tsum(Tsum.group==glab,:);
            if isempty(Ti); continue; end
            % map scene_id -> position index 1..N
            [~, pos] = ismember(double(Ti.scene_id), scenes);
            [pos, o] = sort(pos);
            Ti = Ti(o,:);
            errorbar(pos, Ti.mean, Ti.sem, 'o-', 'LineWidth', 1.8, 'Color', colors(gi,:), ...
                'MarkerFaceColor', colors(gi,:), 'DisplayName', sprintf('%s (mean n≈%.1f)', glab, mean(Ti.n,'omitnan')));
        end

        xlim([0.5, numel(scenes)+0.5]);
        set(gca,'XTick',1:numel(scenes),'XTickLabel',tickLabels);
        xtickangle(30);

        % Block split line between scene 6 and 7
        try
            if numel(scenes) >= 12
                xline(6.5, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, 'HandleVisibility','off');
                yl2 = ylim;
                ytop = yl2(2) - 0.02*(yl2(2)-yl2(1));
                text(3.25, ytop, 'Block 1', 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize', 9);
                text(9.25, ytop, 'Block 2', 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize', 9);
            end
        catch
        end
        ylabel(char(met), 'Interpreter','none');
        xlabel('Scene (design label)', 'Interpreter','none');

        ttag = "";
        if strlength(string(tag))>0
            ttag = " [" + string(tag) + "]";
        end
        title(sprintf('%s: %s across scenes%s', fac, met, ttag), 'Interpreter','none');
        legend('Location','best');

        fn = sprintf('group_scene_sequence_%s_%s', lower(fac), char(met));
        if strlength(string(tag))>0
            fn = fn + "_" + string(tag);
        end
        f1 = fullfile(fp_fig, pipeline.sanitize_filename(fn + ".png"));
        try
            saveas(fig, f1);
        catch
        end
        try; close(fig); catch; end
    end
end

end

function labels = local_scene_labels(T)
labels = repmat("", height(T), 1);
if ismember('scene_name', T.Properties.VariableNames)
    labels = string(T.scene_name);
elseif ismember('SceneID', T.Properties.VariableNames)
    labels = string(T.SceneID);
else
    % fallback: WWR+Cond
    w = repmat("", height(T), 1);
    c = repmat("", height(T), 1);
    if ismember('WWR', T.Properties.VariableNames); w = string(T.WWR); end
    if ismember('Cond', T.Properties.VariableNames); c = string(T.Cond); end
    for i=1:height(T)
        if ~isnan(T.scene_id(i))
            labels(i) = "round" + sprintf('%02d', T.scene_id(i)) + "_W" + w(i) + "_" + c(i);
        end
    end
end
labels = strrep(labels, "_", "-");
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

function out = local_shorten_label(lab)
% Make x-axis labels readable in journal figures.
% Rule: keep alnum + '-' + '_' ; collapse spaces; clip to <=14 chars.
try
    s = string(lab);
    s = strrep(s, " ", "");
    s = strrep(s, "_", "-");
    % remove common prefixes
    s = regexprep(s, '^round\d+[-_]*', '');
    % keep only safe chars
    s = regexprep(s, '[^A-Za-z0-9\-]', '');
    if strlength(s) > 14
        s = extractBetween(s, 1, 14);
    end
    out = s;
catch
    out = string(lab);
end
end
