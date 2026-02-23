function out = compare_denominator_sensitivity(fp_sum, cfg)
%COMPARE_DENOMINATOR_SENSITIVITY Compare 1–30 vs 1–40 relative power results.
%
% This is NOT a full LMM; it mirrors the repository's default aggregation:
% within-subject means by Complexity (Low/High) and between-subject group means.
%
% Outputs:
%   - denominator_sensitivity.csv
%   - denominator_sensitivity_report.md

out = struct();

fp_scene = fullfile(fp_sum, 'all_subjects_scene_level_qc.csv');
if ~exist(fp_scene,'file')
    % fallback to raw
    fp_scene = fullfile(fp_sum, 'all_subjects_scene_level.csv');
end
if ~exist(fp_scene,'file')
    return;
end

S = readtable(fp_scene);
if ~ismember('subject_id', S.Properties.VariableNames) || ~ismember('Complexity', S.Properties.VariableNames)
    return;
end

metrics = {"O_alpha","O_theta","O_beta","F_theta"};
try
    if isfield(cfg,'paper_metrics') && ~isempty(cfg.paper_metrics)
        metrics = cellstr(string(cfg.paper_metrics));
    end
catch
end

% Only metrics with _40 counterpart are comparable
metPairs = {};
for i=1:numel(metrics)
    m = string(metrics{i});
    m40 = m + "_40";
    if ismember(m, S.Properties.VariableNames) && ismember(m40, S.Properties.VariableNames)
        metPairs(end+1,:) = {char(m), char(m40)}; %#ok<AGROW>
    end
end
if isempty(metPairs)
    return;
end

cx01 = local_complexity01(S);
S.cx01 = cx01;

factors = ["Experience","SportFreq"];
Rows = {};

for fi=1:numel(factors)
    fac = factors(fi);
    if ~ismember(fac, S.Properties.VariableNames)
        continue;
    end
    grp = strtrim(string(S.(fac)));
    grp(lower(grp)=="high") = "High";
    grp(lower(grp)=="low") = "Low";

    keep0 = ~isnan(S.cx01) & (grp=="High" | grp=="Low");
    if sum(keep0) < 6
        continue;
    end

    sid = string(S.subject_id);

    for mi=1:size(metPairs,1)
        m30 = string(metPairs{mi,1});
        m40 = string(metPairs{mi,2});

        y30 = double(S.(m30));
        y40 = double(S.(m40));

        keep = keep0 & ~isnan(y30) & ~isnan(y40);
        if sum(keep) < 6
            continue;
        end

        % subject-level mean by Complexity for each denominator
        [G, sid_u, grp_u, cx_u] = findgroups(sid(keep), grp(keep), S.cx01(keep));
        mu30 = splitapply(@(x) mean(x,'omitnan'), y30(keep), G);
        mu40 = splitapply(@(x) mean(x,'omitnan'), y40(keep), G);
        Tsub = table(sid_u, grp_u, cx_u, mu30, mu40, 'VariableNames', {'subject_id','group','cx','mu30','mu40'});

        % deltas per subject
        [Gs, sid2, g2] = findgroups(Tsub.subject_id, Tsub.group);
        d30 = splitapply(@(cx,mu) local_delta(cx,mu), Tsub.cx, Tsub.mu30, Gs);
        d40 = splitapply(@(cx,mu) local_delta(cx,mu), Tsub.cx, Tsub.mu40, Gs);
        D = table(sid2, g2, d30, d40, 'VariableNames', {'subject_id','group','delta30','delta40'});

        % per group summary
        for glab = ["Low","High"]
            Di = D(D.group==glab,:);
            if isempty(Di)
                continue;
            end
            n = sum(~isnan(Di.delta30) & ~isnan(Di.delta40));
            m_d30 = mean(Di.delta30,'omitnan');
            m_d40 = mean(Di.delta40,'omitnan');
            sign_ok = sign(m_d30) == sign(m_d40);
            r = NaN;
            try
                r = corr(Di.delta30, Di.delta40, 'Rows','complete');
            catch
            end

            Rows(end+1,:) = {char(fac), char(glab), char(m30), n, m_d30, m_d40, double(sign_ok), r}; %#ok<AGROW>
        end
    end
end

if isempty(Rows)
    return;
end

T = cell2table(Rows, 'VariableNames', {'factor','group','metric','n_subjects','delta_1_30','delta_1_40','sign_consistent','corr_subject_deltas'});

out.denominator_sensitivity_csv = fullfile(fp_sum, 'denominator_sensitivity.csv');
writetable(T, out.denominator_sensitivity_csv);

% Markdown report
out.denominator_sensitivity_md = fullfile(fp_sum, 'denominator_sensitivity_report.md');
try
    fid = fopen(out.denominator_sensitivity_md,'w');
    if fid~=-1
        fprintf(fid, '# Denominator Sensitivity Report (1–30 vs 1–40)\n\n');
        fprintf(fid, '- Generated: %s\n', datestr(now,31));
        fprintf(fid, '- Note: This report mirrors the repository\'s default aggregation (subject means by Complexity), not a full LMM.\n\n');
        fprintf(fid, 'See `denominator_sensitivity.csv` for the full table.\n');
        fclose(fid);
    end
catch
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

function d = local_delta(cx, mu)
lo = mu(cx==0);
hi = mu(cx==1);
if isempty(lo) || isempty(hi)
    d = NaN;
else
    d = mean(hi,'omitnan') - mean(lo,'omitnan');
end
end
