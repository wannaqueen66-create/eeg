function Stats = compute_blockdiff_stats(S, metrics, mode, gcol)
%COMPUTE_BLOCKDIFF_STATS Helper to compute paired Block2 vs Block1 stats.
% mode:
%   "all"  -> ignore grouping, single row per metric
%   "group"-> split by gcol (Experience/SportFreq) and Group (High/Low)

if nargin < 3 || strlength(string(mode))==0
    mode = "all";
end
if nargin < 4
    gcol = "";
end

metrics = string(metrics);

statsRows = {};

if mode == "all"
    for mi=1:numel(metrics)
        m = metrics(mi);
        X = S(S.metric==m,:);
        if height(X) < 3
            continue;
        end
        y1 = double(X.block1);
        y2 = double(X.block2);
        mask = ~isnan(y1) & ~isnan(y2);
        y1 = y1(mask); y2 = y2(mask);
        n = numel(y1);
        if n < 3
            continue;
        end
        d = y2 - y1;
        md = mean(d,'omitnan');
        sd = std(d,'omitnan');
        sem = sd / sqrt(n);
        dz = md / sd;
        try
            [~, p, ~, st] = ttest(y2, y1);
            tstat = st.tstat; df = st.df;
        catch
            p = NaN; tstat = NaN; df = NaN;
        end
        try
            p_sr = signrank(y2, y1);
        catch
            p_sr = NaN;
        end
        statsRows(end+1,:) = {"All","All", m, n, md, sem, tstat, df, p, dz, p_sr}; %#ok<AGROW>
    end
    Stats = cell2table(statsRows, 'VariableNames', {
        'GroupType','Group','metric','N','mean_diff','sem_diff','t','df','p_ttest','cohen_dz','p_signrank'});
    return;
end

% mode == group
if strlength(string(gcol))==0 || ~ismember(gcol, S.Properties.VariableNames)
    Stats = table();
    return;
end

for g = ["Low","High"]
    for mi=1:numel(metrics)
        m = metrics(mi);
        X = S(S.metric==m & string(S.(gcol))==g, :);
        if height(X) < 3
            continue;
        end
        y1 = double(X.block1);
        y2 = double(X.block2);
        mask = ~isnan(y1) & ~isnan(y2);
        y1 = y1(mask); y2 = y2(mask);
        n = numel(y1);
        if n < 3
            continue;
        end
        d = y2 - y1;
        md = mean(d,'omitnan');
        sd = std(d,'omitnan');
        sem = sd / sqrt(n);
        dz = md / sd;
        try
            [~, p, ~, st] = ttest(y2, y1);
            tstat = st.tstat; df = st.df;
        catch
            p = NaN; tstat = NaN; df = NaN;
        end
        try
            p_sr = signrank(y2, y1);
        catch
            p_sr = NaN;
        end
        statsRows(end+1,:) = {string(gcol), g, m, n, md, sem, tstat, df, p, dz, p_sr}; %#ok<AGROW>
    end
end

Stats = cell2table(statsRows, 'VariableNames', {
    'GroupType','Group','metric','N','mean_diff','sem_diff','t','df','p_ttest','cohen_dz','p_signrank'});
end
