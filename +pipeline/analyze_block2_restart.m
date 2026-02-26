function out = analyze_block2_restart(AllScene, fp_sum, cfg, tag)
%ANALYZE_BLOCK2_RESTART Block2 restart effect: compare Block2 first scene vs mean of remaining 5.
%
% Output root (requested by user): <summary>/analysis-2/
% Writes (per tag: raw/qc):
%   tables/<tag>/block2_restart_subjectlevel_<tag>.csv
%   tables/<tag>/block2_restart_stats_<tag>.csv
%   reports/<tag>/block2_restart_report_<tag>.md
%
% Requires columns in AllScene:
%   subject_id, cond, block_id, cycle_in_block
%   metrics: O_theta, F_theta, O_alpha, O_beta (configurable)
% Optional grouping columns:
%   Experience, SportFreq  (values: High/Low or 高/低 or 1/0)

if nargin < 4 || isempty(tag)
    tag = 'raw';
end
if nargin < 3
    cfg = struct();
end

out = struct();

% --- output dirs ---
fp_a2 = fullfile(fp_sum, 'analysis-2');
fp_tbl = fullfile(fp_a2, 'tables', tag);
fp_rep = fullfile(fp_a2, 'reports', tag);
if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end

% --- required columns ---
req = {'subject_id','cond','block_id','cycle_in_block'};
for i=1:numel(req)
    if ~ismember(req{i}, AllScene.Properties.VariableNames)
        warning('analyze_block2_restart: missing required column %s. Skipping.', req{i});
        return;
    end
end

% metrics list (default requested)
metrics = {"O_theta","F_theta","O_alpha","O_beta"};
try
    if isfield(cfg,'block2_restart_metrics') && ~isempty(cfg.block2_restart_metrics)
        metrics = string(cfg.block2_restart_metrics);
    end
catch
end

% filter Block2 view
T = AllScene;
try
    T.cond = string(T.cond);
end
viewMask = (lower(strtrim(string(T.cond)))=="view");
T = T(viewMask, :);
T = T(double(T.block_id)==2, :);

if isempty(T) || height(T)==0
    warning('analyze_block2_restart: no Block2 view rows found.');
    return;
end

% normalize grouping columns (if present)
if ismember('Experience', T.Properties.VariableNames)
    T.Experience = normalize_high_low(T.Experience);
end
if ismember('SportFreq', T.Properties.VariableNames)
    T.SportFreq = normalize_high_low(T.SportFreq);
end

% subject-level wide to per-metric long
subs = unique(string(T.subject_id), 'stable');

rows = {};
for si=1:numel(subs)
    sid = subs(si);
    Ts = T(string(T.subject_id)==sid, :);

    % First scene in Block2
    iFirst = (double(Ts.cycle_in_block)==1);
    iRest  = (double(Ts.cycle_in_block)>=2 & double(Ts.cycle_in_block)<=6);

    if ~any(iFirst) || ~any(iRest)
        continue;
    end

    % group labels (may be missing)
    ex = ""; sf = "";
    try
        if ismember('Experience', Ts.Properties.VariableNames)
            ex = string(Ts.Experience(find(iFirst,1,'first')));
        end
    catch
    end
    try
        if ismember('SportFreq', Ts.Properties.VariableNames)
            sf = string(Ts.SportFreq(find(iFirst,1,'first')));
        end
    catch
    end

    for mi=1:numel(metrics)
        m = string(metrics(mi));
        if ~ismember(m, Ts.Properties.VariableNames)
            continue;
        end
        yFirst = Ts.(m)(iFirst);
        yRest  = Ts.(m)(iRest);
        yFirst = double(yFirst(1));
        yRestM = mean(double(yRest), 'omitnan');
        if isnan(yFirst) || isnan(yRestM)
            continue;
        end
        diffv = yFirst - yRestM;

        rows(end+1,:) = {sid, m, yFirst, yRestM, diffv, ex, sf}; %#ok<AGROW>
    end
end

if isempty(rows)
    warning('analyze_block2_restart: no usable subject-level pairs found.');
    return;
end

S = cell2table(rows, 'VariableNames', {'subject_id','metric','block2_first','block2_rest_mean','diff','Experience','SportFreq'});

% write subject-level table
fp_sub = fullfile(fp_tbl, sprintf('block2_restart_subjectlevel_%s.csv', tag));
writetable(S, fp_sub);

% stats per groupType x group x metric
statsRows = {};
groupTypes = ["Experience","SportFreq"];
for gi=1:numel(groupTypes)
    gcol = groupTypes(gi);
    if ~ismember(gcol, S.Properties.VariableNames)
        continue;
    end
    for g = ["Low","High"]
        for mi=1:numel(metrics)
            m = string(metrics(mi));
            X = S(S.metric==m & string(S.(gcol))==g, :);
            if height(X) < 3
                continue;
            end
            d = double(X.diff);
            n = sum(~isnan(d));
            if n < 3
                continue;
            end
            md = mean(d,'omitnan');
            sd = std(d,'omitnan');
            sem = sd / sqrt(n);
            dz = md / sd;

            % paired t-test first vs rest
            try
                [~, p, ~, st] = ttest(double(X.block2_first), double(X.block2_rest_mean));
                tstat = st.tstat;
                df = st.df;
            catch
                p = NaN; tstat = NaN; df = NaN;
            end

            % signrank as robust alternative
            try
                p_sr = signrank(double(X.block2_first), double(X.block2_rest_mean));
            catch
                p_sr = NaN;
            end

            statsRows(end+1,:) = {gcol, g, m, n, md, sem, tstat, df, p, dz, p_sr}; %#ok<AGROW>
        end
    end
end

Stats = cell2table(statsRows, 'VariableNames', {
    'GroupType','Group','metric','N','mean_diff','sem_diff','t','df','p_ttest','cohen_dz','p_signrank'});

fp_stats = fullfile(fp_tbl, sprintf('block2_restart_stats_%s.csv', tag));
writetable(Stats, fp_stats);

% report markdown
fp_md = fullfile(fp_rep, sprintf('block2_restart_report_%s.md', tag));
write_report_md(fp_md, tag, metrics, fp_sub, fp_stats, Stats);

out.subjectlevel_csv = fp_sub;
out.stats_csv = fp_stats;
out.report_md = fp_md;
end

function y = normalize_high_low(x)
% Map various encodings -> "High"/"Low" (otherwise empty)
s = lower(strtrim(string(x)));

% common variants
s(ismember(s,["high","h","1","true","yes","y","高"])) = "high";
s(ismember(s,["low","l","0","false","no","n","低"])) = "low";

% final
y = repmat("", numel(s), 1);
y(s=="high") = "High";
y(s=="low")  = "Low";
end

function write_report_md(fp_md, tag, metrics, fp_sub, fp_stats, Stats)
lines = {};
lines{end+1} = sprintf('# Block2 restart analysis (%s)', tag);
lines{end+1} = '';
lines{end+1} = 'Definition:';
lines{end+1} = '- Use **Block2 view** segments only.';
lines{end+1} = '- For each subject and metric: **diff = Block2(first scene) − mean(Block2 scenes 2–6)**.';
lines{end+1} = '';
lines{end+1} = 'Outputs:';
lines{end+1} = sprintf('- Subject-level table: `%s`', fp_sub);
lines{end+1} = sprintf('- Stats table: `%s`', fp_stats);
lines{end+1} = '';
lines{end+1} = 'Metrics:';
lines{end+1} = sprintf('- %s', strjoin(string(metrics), ', '));
lines{end+1} = '';

if isempty(Stats) || height(Stats)==0
    lines{end+1} = '> No group stats computed (missing grouping columns or insufficient N).';
else
    lines{end+1} = 'Quick summary (paired t-test p-values):';
    lines{end+1} = '';
    % write a compact bullet list
    for i=1:height(Stats)
        lines{end+1} = sprintf('- %s %s | %s: N=%d, mean(diff)=%.4g, t(%d)=%.3f, p=%.3g (signrank p=%.3g)', ...
            string(Stats.GroupType(i)), string(Stats.Group(i)), string(Stats.metric(i)), ...
            Stats.N(i), Stats.mean_diff(i), Stats.df(i), Stats.t(i), Stats.p_ttest(i), Stats.p_signrank(i));
    end
end

txt = strjoin(string(lines), newline);
fid = fopen(fp_md,'w');
if fid<0
    warning('analyze_block2_restart: cannot write report %s', fp_md);
    return;
end
fwrite(fid, txt);
fclose(fid);
end
