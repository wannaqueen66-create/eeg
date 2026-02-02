function marker_interval_stats(inputPath, varargin)
%MARKER_INTERVAL_STATS  Export consecutive marker intervals for manual QC/fixing.
%
% Usage:
%   marker_interval_stats('/path/to/file.set')
%   marker_interval_stats('/path/to/folder_with_set_files')
%
% Output (per dataset):
%   <base>_marker_intervals.csv              (row per consecutive marker pair)
%   <base>_marker_transition_summary.csv     (stats grouped by m0->m1)
%
% Notes:
% - Assumes EEGLAB is available (pop_loadset). Add EEGLAB to MATLAB path first.
% - Markers are expected to be numeric 1-9, stored in EEG.event.type.
%
% Optional args (name-value):
%   'OutDir'        : output directory (default: alongside dataset)
%   'MarkerRange'   : vector of allowed markers (default: 1:9)
%   'Recursive'     : when inputPath is folder, search recursively (default false)

p = inputParser;
p.addRequired('inputPath', @(x) ischar(x) || isstring(x));
p.addParameter('OutDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('MarkerRange', 1:9, @(x) isnumeric(x) && isvector(x));
p.addParameter('Recursive', false, @(x) islogical(x) && isscalar(x));
p.parse(inputPath, varargin{:});
args = p.Results;

inputPath = char(inputPath);
markerRange = args.MarkerRange;

if isfolder(inputPath)
    % collect .set files
    if args.Recursive
        files = dir(fullfile(inputPath, '**', '*.set'));
    else
        files = dir(fullfile(inputPath, '*.set'));
    end
    assert(~isempty(files), 'No .set files found in folder: %s', inputPath);

    % aggregate transition stats across all files
    agg = []; % struct array for per-transition rows

    for i = 1:numel(files)
        fp = fullfile(files(i).folder, files(i).name);
        fprintf('[%d/%d] %s\n', i, numel(files), fp);
        out = process_one_set(fp, args.OutDir, markerRange);
        if ~isempty(out.summaryRows)
            agg = [agg; out.summaryRows]; %#ok<AGROW>
        end
    end

    % write aggregated summary
    if ~isempty(agg)
        T = struct2table(agg);
        % group again across all files
        [G, m0, m1] = findgroups(T.m0, T.m1);
        cnt = splitapply(@sum, T.count, G);
        % for mean/median across files we recompute using expanded stats is hard; instead
        % we approximate by weighting file-level mean by count.
        wmean = splitapply(@(x,w) sum(x.*w)./sum(w), T.mean_s, T.count, G);
        wstd  = splitapply(@(x,w) sqrt(sum(((x - sum(x.*w)./sum(w)).^2).*w)./sum(w)), T.mean_s, T.count, G);
        minv  = splitapply(@min, T.min_s, G);
        maxv  = splitapply(@max, T.max_s, G);

        Tout = table(m0, m1, cnt, wmean, wstd, minv, maxv, ...
            'VariableNames', {'m0','m1','count','wmean_s','wstd_s','min_s','max_s'});

        if strlength(string(args.OutDir)) > 0
            outdir = char(args.OutDir);
        else
            outdir = inputPath;
        end
        if ~exist(outdir, 'dir'); mkdir(outdir); end
        fp_out = fullfile(outdir, 'ALL_marker_transition_summary.csv');
        writetable(Tout, fp_out);
        fprintf('[OK] Wrote aggregated transition summary: %s\n', fp_out);
    end

else
    process_one_set(inputPath, args.OutDir, markerRange);
end

end

function out = process_one_set(fp_set, outDir, markerRange)
% Load dataset
assert(exist('pop_loadset', 'file') == 2, 'EEGLAB function pop_loadset not found. Add EEGLAB to MATLAB path first.');
EEG = pop_loadset('filename', fp_set);
fs = double(EEG.srate);

% extract markers from EEG.event
if ~isfield(EEG, 'event') || isempty(EEG.event)
    warning('No EEG.event in dataset: %s', fp_set);
    out.summaryRows = [];
    return;
end

types_raw = {EEG.event.type};
lat_raw   = [EEG.event.latency];

mk = nan(size(types_raw));
for i = 1:numel(types_raw)
    if isnumeric(types_raw{i})
        mk(i) = double(types_raw{i});
    else
        v = str2double(string(types_raw{i}));
        if ~isnan(v), mk(i) = v; end
    end
end

keep = ismember(mk, markerRange);
mk   = mk(keep);
lat  = double(lat_raw(keep));

% sort by latency
[lat, ord] = sort(lat);
mk = mk(ord);

n = numel(mk);
if n < 2
    warning('Not enough valid markers (%d) in dataset: %s', n, fp_set);
    out.summaryRows = [];
    return;
end

% consecutive intervals
m0 = mk(1:end-1)';
m1 = mk(2:end)';
lat0 = lat(1:end-1)';
lat1 = lat(2:end)';

t0 = (lat0-1) / fs;
t1 = (lat1-1) / fs;
dt = t1 - t0;

idx = (1:numel(dt))';
T = table(idx, m0, m1, lat0, lat1, t0, t1, dt, ...
    'VariableNames', {'idx','m0','m1','lat0_samp','lat1_samp','t0_s','t1_s','dt_s'});

% transition summary
[G, gm0, gm1] = findgroups(T.m0, T.m1);
count = splitapply(@numel, T.dt_s, G);
mean_s = splitapply(@mean, T.dt_s, G);
median_s = splitapply(@median, T.dt_s, G);
std_s = splitapply(@std, T.dt_s, G);
min_s = splitapply(@min, T.dt_s, G);
max_s = splitapply(@max, T.dt_s, G);
Tsum = table(gm0, gm1, count, mean_s, median_s, std_s, min_s, max_s, ...
    'VariableNames', {'m0','m1','count','mean_s','median_s','std_s','min_s','max_s'});

% output paths
[folder, base, ~] = fileparts(fp_set);
if strlength(string(outDir)) > 0
    od = char(outDir);
else
    od = folder;
end
if ~exist(od, 'dir'); mkdir(od); end

fp_intervals = fullfile(od, sprintf('%s_marker_intervals.csv', base));
fp_summary   = fullfile(od, sprintf('%s_marker_transition_summary.csv', base));

writetable(T, fp_intervals);
writetable(Tsum, fp_summary);

fprintf('[OK] Wrote: %s\n', fp_intervals);
fprintf('[OK] Wrote: %s\n', fp_summary);

% return rows for aggregation
summaryRows = table2struct(Tsum);
for k = 1:numel(summaryRows)
    summaryRows(k).dataset = string(base);
end
out.summaryRows = summaryRows;

end
