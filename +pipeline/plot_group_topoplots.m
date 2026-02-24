function plot_group_topoplots(fp_out, fp_sum, cfg)
%PLOT_GROUP_TOPOPLOTS Group-level topoplots from per-subject exported topo CSV.
%
% Requires per-subject files:
%   <fp_out>/<subject>/csv/<subject>_topo_long.csv
% and at least one chanlocs snapshot:
%   <fp_out>/<subject>/qc/<subject>_chanlocs.mat
%
% Produces PNGs under:
%   <fp_sum>/fig/
%
% Factors (between-subject): Experience, SportFreq (High/Low)
% Metric (within-subject): viewComplexityHigh_minus_viewComplexityLow
% Bands: theta, alpha, beta

if nargin < 3
    cfg = struct();
end

% Require EEGLAB topoplot
if exist('topoplot','file') ~= 2
    warning('plot_group_topoplots: topoplot not found on MATLAB path. Start EEGLAB or add it to path, then rerun summarize_bandpower_outputs.');
    return;
end

tag = 'raw';
try
    if isfield(cfg,'qc_include_subjects') && ~isempty(cfg.qc_include_subjects)
        tag = 'qc';
    end
catch
end
fp_fig = pipeline.get_fig_dir(fp_sum, cfg, 'topo', tag);

% Load chanlocs from first available snapshot
D = dir(fullfile(fp_out, '*', 'qc', '*_chanlocs.mat'));
if isempty(D)
    warning('plot_group_topoplots: no *_chanlocs.mat found. Skipping group topoplots.');
    return;
end
S = load(fullfile(D(1).folder, D(1).name));
if ~isfield(S,'chanlocs')
    warning('plot_group_topoplots: chanlocs variable not found in %s', fullfile(D(1).folder, D(1).name));
    return;
end
chanlocs = S.chanlocs;
nbchan = numel(chanlocs);

% Read all per-subject topo exports
Tall = table();
subs = dir(fp_out);
subs = subs([subs.isdir]);
subs = setdiff({subs.name},{'.','..','summary'});
subs = sort(subs);

for k=1:numel(subs)
    sid = string(subs{k});
    f = fullfile(fp_out, char(sid), 'csv', sprintf('%s_topo_long.csv', sid));
    if exist(f,'file')
        try
            T = readtable(f, 'TextType','string');
            Tall = [Tall; T]; %#ok<AGROW>
        catch ME
            fprintf(2, '[WARN] Failed reading topo_long for %s: %s\n', sid, ME.message);
        end
    end
end

if isempty(Tall) || height(Tall)==0
    warning('plot_group_topoplots: no topo_long data found.');
    return;
end

% Apply QC inclusion list (optional)
try
    if isfield(cfg,'qc_include_subjects') && ~isempty(cfg.qc_include_subjects)
        inc = string(cfg.qc_include_subjects);
        Tall = Tall(ismember(string(Tall.subject_id), inc), :);
    end
catch
end

% Ensure string columns for robust comparisons (case-insensitive matching)
try
    if ismember('metric', Tall.Properties.VariableNames)
        Tall.metric = strtrim(string(Tall.metric));
    end
    if ismember('band', Tall.Properties.VariableNames)
        Tall.band = lower(strtrim(string(Tall.band)));
    end
catch
end

% Attach between-subject factors (Experience/SportFreq) if missing in topo_long
% topo_long exported from single-subject stage may only include subject_id/chan/band/metric/value.
if ~ismember('Experience', Tall.Properties.VariableNames) || ~ismember('SportFreq', Tall.Properties.VariableNames)
    try
        map = table();

        % Tidy layout: tables live under summary/tables/merged_raw
        fp_tbl_raw = fp_sum;
        try
            if exist('pipeline.get_table_dir','file')==2
                fp_tbl_raw = pipeline.get_table_dir(fp_sum, cfg, 'merged_raw');
            end
        catch
        end

        f_map = fullfile(fp_tbl_raw, 'per_subject_recovery_metrics.csv');
        if exist(f_map,'file')
            map = readtable(f_map, 'TextType','string');
        else
            f_pairs = fullfile(fp_tbl_raw, 'all_subjects_pairs_check.csv');
            if exist(f_pairs,'file')
                map = readtable(f_pairs, 'TextType','string');
            else
                % legacy fallback
                f_map2 = fullfile(fp_sum, 'per_subject_recovery_metrics.csv');
                if exist(f_map2,'file')
                    map = readtable(f_map2, 'TextType','string');
                else
                    f_pairs2 = fullfile(fp_sum, 'all_subjects_pairs_check.csv');
                    if exist(f_pairs2,'file')
                        map = readtable(f_pairs2, 'TextType','string');
                    end
                end
            end
        end

        if ~isempty(map) && ismember('subject_id', map.Properties.VariableNames)
            % Reduce to unique subject rows
            sidm = string(map.subject_id);
            [~, ia] = unique(sidm, 'stable');
            map = map(ia,:);

            if ~ismember('SportFreq', Tall.Properties.VariableNames)
                sf = repmat("", height(Tall), 1);
                if ismember('SportFreq', map.Properties.VariableNames)
                    [tf, loc] = ismember(string(Tall.subject_id), string(map.subject_id));
                    sf(tf) = strtrim(string(map.SportFreq(loc(tf))));
                end
                Tall.SportFreq = sf;
            end

            if ~ismember('Experience', Tall.Properties.VariableNames)
                ex = repmat("", height(Tall), 1);
                if ismember('Experience', map.Properties.VariableNames)
                    [tf, loc] = ismember(string(Tall.subject_id), string(map.subject_id));
                    ex(tf) = strtrim(string(map.Experience(loc(tf))));
                end
                Tall.Experience = ex;
            end
        end
    catch
    end
end

% Normalize factor columns
if ismember('Experience', Tall.Properties.VariableNames)
    Tall.Experience = strtrim(string(Tall.Experience));
end
if ismember('SportFreq', Tall.Properties.VariableNames)
    Tall.SportFreq = strtrim(string(Tall.SportFreq));
end

% Fallback: if factors still missing/empty, try mapping from all_subjects_scene_level
try
    needEx = ~ismember('Experience', Tall.Properties.VariableNames) || all(strlength(strtrim(string(Tall.Experience)))==0);
    needSf = ~ismember('SportFreq', Tall.Properties.VariableNames) || all(strlength(strtrim(string(Tall.SportFreq)))==0);
    if needEx || needSf
        fp_tbl_raw = fp_sum;
        try
            if exist('pipeline.get_table_dir','file')==2
                fp_tbl_raw = pipeline.get_table_dir(fp_sum, cfg, 'merged_raw');
            end
        catch
        end

        f_scene = fullfile(fp_tbl_raw, 'all_subjects_scene_level.csv');
        if ~exist(f_scene,'file')
            % legacy fallback
            f_scene = fullfile(fp_sum, 'all_subjects_scene_level.csv');
        end

        if exist(f_scene,'file')
            M = readtable(f_scene, 'TextType','string');
            if ismember('subject_id', M.Properties.VariableNames)
                sidm = string(M.subject_id);
                [~, ia] = unique(sidm, 'stable');
                M = M(ia,:);

                if needEx && ismember('Experience', M.Properties.VariableNames)
                    [tf, loc] = ismember(string(Tall.subject_id), string(M.subject_id));
                    ex = repmat("", height(Tall), 1);
                    ex(tf) = strtrim(string(M.Experience(loc(tf))));
                    Tall.Experience = ex;
                end
                if needSf && ismember('SportFreq', M.Properties.VariableNames)
                    [tf, loc] = ismember(string(Tall.subject_id), string(M.subject_id));
                    sf = repmat("", height(Tall), 1);
                    sf(tf) = strtrim(string(M.SportFreq(loc(tf))));
                    Tall.SportFreq = sf;
                end
            end
        end
    end
catch
end

% If still missing factor columns, bail with a clear warning (otherwise silently produces no plots)
if ~ismember('Experience', Tall.Properties.VariableNames) && ~ismember('SportFreq', Tall.Properties.VariableNames)
    warning('plot_group_topoplots: Experience/SportFreq not found in topo_long and could not be attached from summary tables. No group topoplots were generated.');
    return;
end

% Defaults
metricName = "viewComplexityHigh_minus_viewComplexityLow";
bands = ["theta","alpha","beta"];

factors = ["Experience","SportFreq"];

for fi=1:numel(factors)
    fac = factors(fi);
    if ~ismember(fac, Tall.Properties.VariableNames)
        continue;
    end

    for bi=1:numel(bands)
        band = bands(bi);

        m = strtrim(string(Tall.metric));
        % metric match (case-insensitive)
        use = (lower(m)==lower(metricName));
        % band match (already lowercased)
        use = use & (Tall.band==band);
        % normalize factor values to High/Low (support Chinese/1-0)
        fval = lower(strtrim(string(Tall.(fac))));
        fval(fval=="高" | fval=="1" | fval=="h") = "high";
        fval(fval=="低" | fval=="0" | fval=="l") = "low";
        use = use & (fval=="high" | fval=="low");

        if sum(use) < nbchan*4
            fprintf(2,'[WARN] topo: not enough rows for %s %s (%s). got=%d need>=%d\n', fac, band, metricName, sum(use), nbchan*4);
            continue;
        end

        X = Tall(use,:);

        % Ensure required columns
        req = {'subject_id','chan_idx','value'};
        for r=1:numel(req)
            if ~ismember(req{r}, X.Properties.VariableNames)
                warning('plot_group_topoplots: missing column %s in topo_long.csv', req{r});
                return;
            end
        end

        % Standardize group labels
        grp = strtrim(string(X.(fac)));
        grp(lower(grp)=="high") = "High";
        grp(lower(grp)=="low") = "Low";
        X.group = grp;

        % Compute subject-level vectors first (robust)
        % One value per (subject, chan)
        [G1, sid_u, grp_u, ch_u] = findgroups(string(X.subject_id), X.group, double(X.chan_idx));
        mu = splitapply(@(v) mean(double(v),'omitnan'), X.value, G1);
        Ssub = table(sid_u, grp_u, ch_u, mu, 'VariableNames', {'subject_id','group','chan_idx','mu'});

        % Group mean per channel
        groups = ["Low","High"];
        topo = struct();
        for gi=1:numel(groups)
            gname = groups(gi);
            G = Ssub(Ssub.group==gname,:);
            if isempty(G)
                continue;
            end
            vec = nan(nbchan,1);
            for ch=1:nbchan
                vv = G.mu(G.chan_idx==ch);
                if ~isempty(vv)
                    vec(ch) = mean(vv,'omitnan');
                end
            end
            topo.(char(gname)) = vec;
        end

                if ~isfield(topo,'Low') || ~isfield(topo,'High')
            fprintf(2,'[WARN] topo: missing Low/High for %s %s after grouping\n', fac, band);
            continue;
        end

        % Plot Low
        fig = figure('Name', sprintf('%s %s %s Low', fac, band, metricName), 'Position',[100 100 700 520]);
        topoplot(topo.Low, chanlocs, 'electrodes','labels');
        title(sprintf('%s | %s | %s | Low', fac, band, metricName), 'Interpreter','none');
        colorbar;
        f1 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('group_topo_%s_%s_%s_Low.png', lower(fac), band, metricName)));
        saveas(fig, f1);
        try; close(fig); catch; end

        % Plot High
        fig = figure('Name', sprintf('%s %s %s High', fac, band, metricName), 'Position',[100 100 700 520]);
        topoplot(topo.High, chanlocs, 'electrodes','labels');
        title(sprintf('%s | %s | %s | High', fac, band, metricName), 'Interpreter','none');
        colorbar;
        f2 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('group_topo_%s_%s_%s_High.png', lower(fac), band, metricName)));
        saveas(fig, f2);
        try; close(fig); catch; end

        % Plot High - Low
        fig = figure('Name', sprintf('%s %s %s High-Low', fac, band, metricName), 'Position',[100 100 700 520]);
        topoplot(topo.High - topo.Low, chanlocs, 'electrodes','labels');
        title(sprintf('%s | %s | %s | High - Low', fac, band, metricName), 'Interpreter','none');
        colorbar;
        f3 = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('group_topo_%s_%s_%s_HighMinusLow.png', lower(fac), band, metricName)));
        saveas(fig, f3);
        try; close(fig); catch; end

    end
end

end
