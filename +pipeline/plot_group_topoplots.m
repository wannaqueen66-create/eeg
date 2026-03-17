function plot_group_topoplots(fp_out, fp_batch, cfg)
%PLOT_GROUP_TOPOPLOTS Group-level topoplots from per-subject exported topo CSV.
%
% Requires per-subject files:
%   <fp_out>/<subject>/csv/<subject>_topo_long.csv
% and at least one chanlocs snapshot:
%   <fp_out>/<subject>/qc/<subject>_chanlocs.mat
%
% Produces PNGs under:
%   <fp_batch>/fig/
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
fp_fig = pipeline.get_fig_dir(fp_batch, cfg, 'topo', tag);

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

% Canonical subject id (robust join key)
try
    Tall.sid_key = canonical_subject_id(Tall.subject_id);
catch
    Tall.sid_key = strtrim(string(Tall.subject_id));
end

% Apply QC inclusion list (optional)
try
    if isfield(cfg,'qc_include_subjects') && ~isempty(cfg.qc_include_subjects)
        inc = canonical_subject_id(cfg.qc_include_subjects);
        Tall = Tall(ismember(Tall.sid_key, inc), :);
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

% Attach between-subject factors (ExperienceGroup preferred) if missing in topo_long
% topo_long exported from single-subject stage may only include subject_id/chan/band/metric/value.
% First priority: direct group map produced in current summarize run.
if (~ismember('ExperienceGroup', Tall.Properties.VariableNames) || all(strlength(strtrim(string(Tall.ExperienceGroup)))==0)) || ...
   (~ismember('SportFreqGroup', Tall.Properties.VariableNames) || all(strlength(strtrim(string(Tall.SportFreqGroup)))==0))
    try
        if isfield(cfg,'topo_group_map_path') && ~isempty(cfg.topo_group_map_path) && exist(cfg.topo_group_map_path,'file')
            M = readtable(cfg.topo_group_map_path, 'TextType','string');
            if ismember('subject_id', M.Properties.VariableNames)
                sidm = canonical_subject_id(M.subject_id);
                [~, ia] = unique(sidm,'stable');
                M = M(ia,:);
                sidm = sidm(ia);
                [tf, loc] = ismember(Tall.sid_key, sidm);

                if ~ismember('ExperienceGroup', Tall.Properties.VariableNames)
                    Tall.ExperienceGroup = repmat("", height(Tall), 1);
                end
                if ~ismember('SportFreqGroup', Tall.Properties.VariableNames)
                    Tall.SportFreqGroup = repmat("", height(Tall), 1);
                end

                if ismember('ExperienceGroup', M.Properties.VariableNames)
                    tmp = Tall.ExperienceGroup;
                    tmp(tf) = normalize_high_low_local(M.ExperienceGroup(loc(tf)));
                    Tall.ExperienceGroup = tmp;
                end
                if ismember('SportFreqGroup', M.Properties.VariableNames)
                    tmp = Tall.SportFreqGroup;
                    tmp(tf) = normalize_high_low_local(M.SportFreqGroup(loc(tf)));
                    Tall.SportFreqGroup = tmp;
                end
            end
        end
    catch
    end
end

if ~ismember('Experience', Tall.Properties.VariableNames) || ~ismember('SportFreq', Tall.Properties.VariableNames)
    try
        map = table();

        % Tidy layout: tables live under batch/tables/merged_raw
        fp_tbl_raw = fp_batch;
        try
            if exist('pipeline.get_table_dir','file')==2
                fp_tbl_raw = pipeline.get_table_dir(fp_batch, cfg, 'merged_raw');
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
                f_map2 = fullfile(fp_batch, 'per_subject_recovery_metrics.csv');
                if exist(f_map2,'file')
                    map = readtable(f_map2, 'TextType','string');
                else
                    f_pairs2 = fullfile(fp_batch, 'all_subjects_pairs_check.csv');
                    if exist(f_pairs2,'file')
                        map = readtable(f_pairs2, 'TextType','string');
                    end
                end
            end
        end

        if ~isempty(map) && ismember('subject_id', map.Properties.VariableNames)
            sidm = canonical_subject_id(map.subject_id);
            [~, ia] = unique(sidm, 'stable');
            map = map(ia,:);
            sidm = sidm(ia);

            if ~ismember('SportFreqGroup', Tall.Properties.VariableNames)
                sf = repmat("", height(Tall), 1);
                if ismember('SportFreqGroup', map.Properties.VariableNames)
                    [tf, loc] = ismember(Tall.sid_key, sidm);
                    sf(tf) = strtrim(string(map.SportFreqGroup(loc(tf))));
                end
                Tall.SportFreqGroup = sf;
            end

            if ~ismember('ExperienceGroup', Tall.Properties.VariableNames)
                ex = repmat("", height(Tall), 1);
                if ismember('ExperienceGroup', map.Properties.VariableNames)
                    [tf, loc] = ismember(Tall.sid_key, sidm);
                    ex(tf) = strtrim(string(map.ExperienceGroup(loc(tf))));
                end
                Tall.ExperienceGroup = ex;
            end
        end
    catch
    end
end

% Canonical group columns
if ismember('ExperienceGroup', Tall.Properties.VariableNames)
    Tall.ExperienceGroup = strtrim(string(Tall.ExperienceGroup));
end
if ismember('SportFreqGroup', Tall.Properties.VariableNames)
    Tall.SportFreqGroup = strtrim(string(Tall.SportFreqGroup));
end

% Fallback: if canonical group columns still missing/empty, map from all_subjects_scene_level
try
    needExG = ~ismember('ExperienceGroup', Tall.Properties.VariableNames) || all(strlength(strtrim(string(Tall.ExperienceGroup)))==0);
    needSfG = ~ismember('SportFreqGroup', Tall.Properties.VariableNames) || all(strlength(strtrim(string(Tall.SportFreqGroup)))==0);
    if needExG || needSfG
        fp_tbl_raw = fp_batch;
        try
            if exist('pipeline.get_table_dir','file')==2
                fp_tbl_raw = pipeline.get_table_dir(fp_batch, cfg, 'merged_raw');
            end
        catch
        end

        f_scene = fullfile(fp_tbl_raw, 'all_subjects_scene_level.csv');
        if ~exist(f_scene,'file')
            % legacy fallback
            f_scene = fullfile(fp_batch, 'all_subjects_scene_level.csv');
        end

        if exist(f_scene,'file')
            M = readtable(f_scene, 'TextType','string');
            if ismember('subject_id', M.Properties.VariableNames)
                sidm = canonical_subject_id(M.subject_id);
                [~, ia] = unique(sidm, 'stable');
                M = M(ia,:);
                sidm = sidm(ia);

                [tf, loc] = ismember(Tall.sid_key, sidm);

                if needExG
                    ex = repmat("", height(Tall), 1);
                    src = string([]);
                    if ismember('ExperienceGroup', M.Properties.VariableNames)
                        src = string(M.ExperienceGroup(loc(tf)));
                    elseif ismember('Experience', M.Properties.VariableNames)
                        src = string(M.Experience(loc(tf)));
                    end
                    if ~isempty(src)
                        ex(tf) = normalize_high_low_local(src);
                    end
                    Tall.ExperienceGroup = ex;
                end

                if needSfG
                    sf = repmat("", height(Tall), 1);
                    src = string([]);
                    if ismember('SportFreqGroup', M.Properties.VariableNames)
                        src = string(M.SportFreqGroup(loc(tf)));
                    elseif ismember('SportFreq', M.Properties.VariableNames)
                        src = string(M.SportFreq(loc(tf)));
                    end
                    if ~isempty(src)
                        sf(tf) = normalize_high_low_local(src);
                    end
                    Tall.SportFreqGroup = sf;
                end
            end
        end
    end
catch
end

% If still missing or invalid factor columns, bail with a clear warning
hasEx = ismember('ExperienceGroup', Tall.Properties.VariableNames);
validEx = false;
if hasEx
    ex = lower(strtrim(string(Tall.ExperienceGroup)));
    validEx = any(ex=="high" | ex=="low" | ex=="高" | ex=="低" | ex=="1" | ex=="0" | ex=="h" | ex=="l");
end
if ~validEx
    try
        nSubTopo = numel(unique(Tall.sid_key));
    catch
        nSubTopo = NaN;
    end
    warning('plot_group_topoplots: ExperienceGroup not found in topo_long and could not be attached from summary tables. No group topoplots were generated. topo_subjects=%g', nSubTopo);
    return;
end

% Defaults
metricName = "viewComplexityHigh_minus_viewComplexityLow";
bands = ["theta","alpha","beta"];

factors = ["ExperienceGroup"];

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

function g = normalize_high_low_local(x)
s = lower(strtrim(string(x)));
g = repmat("", numel(s), 1);

% strict matches
g(s=="high" | s=="1" | s=="高" | s=="h") = "High";
g(s=="low"  | s=="0" | s=="低" | s=="l") = "Low";

% relaxed contains (e.g., 高组/低组/high group/low group)
maskEmpty = (g=="");
g(maskEmpty & (contains(s,"high") | contains(s,"高"))) = "High";
g(maskEmpty & (contains(s,"low")  | contains(s,"低"))) = "Low";

t = strtrim(string(x));
g(t=="High") = "High";
g(t=="Low") = "Low";
end

function sid = canonical_subject_id(x)
sid = strtrim(string(x));
sid = replace(sid, "\", "/");
% keep basename if accidental path appears
for i=1:numel(sid)
    parts = split(sid(i), "/");
    sid(i) = parts(end);
end
sid = regexprep(sid, '\\.set$', '', 'ignorecase');
sid = strtrim(sid);
end
