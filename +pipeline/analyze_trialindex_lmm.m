function out = analyze_trialindex_lmm(AllScene, fp_sum, cfg, tag)
%ANALYZE_TRIALINDEX_LMM Fit LMM with TrialIndex as covariate (adaptation trend).
%
% Model (per metric, per grouping type):
%   EEG ~ Group*WWR*Complexity + TrialIndex + Group:TrialIndex + (1|Subject)
%
% Uses view segments only.
%
% Outputs under:
%   <summary>/analysis-2/task3_trialindex_lmm/
%     tables/<tag>/{experience|sportfreq}/analysis_ready_<metric>_<tag>.csv
%     tables/<tag>/{experience|sportfreq}/lmm_fixed_effects_<metric>_<tag>.csv
%     tables/<tag>/{experience|sportfreq}/lmm_anova_<metric>_<tag>.csv
%     reports/<tag>/{experience|sportfreq}/trialindex_lmm_report_<metric>_<tag>.md
%
% Notes:
% - Requires Statistics and Machine Learning Toolbox (fitlme). If unavailable,
%   writes a warning and skips.
% - Requires design columns: WWR and Complexity, plus group columns.

if nargin < 4 || isempty(tag)
    tag = 'raw';
end
if nargin < 3
    cfg = struct();
end

out = struct();

% Output dirs
fp_root = fullfile(fp_sum, 'analysis-2', 'task3_trialindex_lmm');
fp_tbl = fullfile(fp_root, 'tables', tag);
fp_rep = fullfile(fp_root, 'reports', tag);
if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end

% Check fitlme
if exist('fitlme','file') ~= 2
    warning('analyze_trialindex_lmm: fitlme not found (Stats toolbox missing). Skipping.');
    return;
end

% Required columns
req = {'subject_id','scene_id'};
for i=1:numel(req)
    if ~ismember(req{i}, AllScene.Properties.VariableNames)
        warning('analyze_trialindex_lmm: missing required column %s. Skipping.', req{i});
        return;
    end
end

% Metrics
metrics = {"O_theta","F_theta","O_alpha","O_beta"};
try
    if isfield(cfg,'trialindex_lmm_metrics') && ~isempty(cfg.trialindex_lmm_metrics)
        metrics = string(cfg.trialindex_lmm_metrics);
    end
catch
end

T = AllScene;
% scene_level tables are view-only; `cond` may not exist here
try
    if ismember('cond', T.Properties.VariableNames)
        T.cond = string(T.cond);
        T = T(lower(strtrim(string(T.cond)))=="view", :);
    end
catch
end

% TrialIndex: default use scene_id (1..12)
TrialIndex = double(T.scene_id);
% if missing/invalid, try from block_id + cycle_in_block
if all(isnan(TrialIndex)) || any(TrialIndex==0)
    if ismember('block_id', T.Properties.VariableNames) && ismember('cycle_in_block', T.Properties.VariableNames)
        TrialIndex = (double(T.block_id)-1)*6 + double(T.cycle_in_block);
    end
end
T.TrialIndex = TrialIndex;

% Normalize WWR / Complexity if present
if ismember('WWR', T.Properties.VariableNames)
    T.WWR = normalize_wwr(T.WWR);
end
if ismember('Complexity', T.Properties.VariableNames)
    T.Complexity = normalize_complexity(T.Complexity);
end

% Group columns normalization (group labels are the canonical grouping vars)
if ismember('ExperienceGroup', T.Properties.VariableNames)
    T.ExperienceGroup = normalize_high_low(T.ExperienceGroup);
end
if ismember('SportFreqGroup', T.Properties.VariableNames)
    T.SportFreqGroup = normalize_high_low(T.SportFreqGroup);
end

% Remove rows with missing essentials
baseNeed = ~isnan(T.TrialIndex) & T.TrialIndex>=1 & T.TrialIndex<=12;
T = T(baseNeed,:);

if height(T) == 0
    warning('analyze_trialindex_lmm: no usable view rows with TrialIndex in [1,12].');
    return;
end

% Setup grouping analyses
analyses = { ...
    struct('name','experience','gcol','ExperienceGroup'), ...
    struct('name','sportfreq','gcol','SportFreqGroup') ...
};

for ai=1:numel(analyses)
    A = analyses{ai};
    gcol = string(A.gcol);
    if ~ismember(gcol, T.Properties.VariableNames)
        continue;
    end

    % Need WWR + Complexity for the requested model
    if ~ismember('WWR', T.Properties.VariableNames) || ~ismember('Complexity', T.Properties.VariableNames)
        warning('analyze_trialindex_lmm: missing WWR/Complexity columns. Provide design_path in run_eeg_bandpower_pipeline so they can be attached.');
        return;
    end

    % Keep only High/Low
    grp = string(T.(gcol));
    use = (grp=="High" | grp=="Low");
    Ta = T(use,:);

    if height(Ta) < 20
        continue;
    end

    % Categorical encoding
    Ta.Subject = categorical(string(Ta.subject_id));
    Ta.Group = categorical(string(Ta.(gcol)));
    Ta.WWR = categorical(string(Ta.WWR));
    Ta.Complexity = categorical(string(Ta.Complexity));

    % prepare dirs
    fp_tbl2 = fullfile(fp_tbl, char(A.name));
    fp_rep2 = fullfile(fp_rep, char(A.name));
    fp_fig2 = fullfile(fp_root, 'figures', tag, char(A.name));
    if ~exist(fp_tbl2,'dir'); mkdir(fp_tbl2); end
    if ~exist(fp_rep2,'dir'); mkdir(fp_rep2); end
    if ~exist(fp_fig2,'dir'); mkdir(fp_fig2); end

    for mi=1:numel(metrics)
        m = string(metrics(mi));
        if ~ismember(m, Ta.Properties.VariableNames)
            continue;
        end

        Y = double(Ta.(m));
        ok = ~isnan(Y);
        Tm = Ta(ok,:);
        Tm.Y = double(Tm.(m));

        if height(Tm) < 30
            continue;
        end

        % analysis-ready long table
        Tready = table(string(Tm.subject_id), Tm.Subject, Tm.Group, Tm.WWR, Tm.Complexity, double(Tm.TrialIndex), double(Tm.Y), ...
            'VariableNames', {'subject_id','Subject','Group','WWR','Complexity','TrialIndex','EEG'});
        fp_ready = fullfile(fp_tbl2, pipeline.sanitize_filename(sprintf('analysis_ready_%s_%s.csv', m, tag)));
        writetable(Tready, fp_ready);

        % Fit LMM
        formula = 'EEG ~ Group*WWR*Complexity + TrialIndex + Group:TrialIndex + (1|Subject)';
        try
            lme = fitlme(Tready, formula);
        catch ME
            fprintf(2, '[WARN] analyze_trialindex_lmm: fitlme failed for %s (%s): %s\n', m, A.name, ME.message);
            continue;
        end

        % Fixed effects table
        try
            FE = lme.Coefficients;
            fp_fe = fullfile(fp_tbl2, pipeline.sanitize_filename(sprintf('lmm_fixed_effects_%s_%s.csv', m, tag)));
            writetable(FE, fp_fe);
        catch
            fp_fe = "";
        end

        % ANOVA table
        try
            AN = anova(lme);
            fp_an = fullfile(fp_tbl2, pipeline.sanitize_filename(sprintf('lmm_anova_%s_%s.csv', m, tag)));
            writetable(AN, fp_an);
        catch
            AN = table();
            fp_an = "";
        end

        % Paper-friendly figures: trend over TrialIndex by Group (overall + faceted by WWR x Complexity)
        try
            pipeline.plot_trialindex_trend(Tready, A.name, m, tag, fp_fig2, lme);
        catch ME
            fprintf(2, '[WARN] analyze_trialindex_lmm: plot_trialindex_trend failed for %s (%s): %s\n', m, A.name, ME.message);
        end
        try
            pipeline.plot_trialindex_trend_facets(Tready, A.name, m, tag, fp_fig2, lme);
        catch ME
            fprintf(2, '[WARN] analyze_trialindex_lmm: plot_trialindex_trend_facets failed for %s (%s): %s\n', m, A.name, ME.message);
        end

        % Report md (highlight TrialIndex + Group:TrialIndex)
        fp_md = fullfile(fp_rep2, pipeline.sanitize_filename(sprintf('trialindex_lmm_report_%s_%s.md', m, tag)));
        write_report(fp_md, A.name, m, tag, fp_ready, fp_fe, fp_an, lme, AN, fp_fig2);

        out.(char(A.name)).(char(m)).ready = fp_ready;
        out.(char(A.name)).(char(m)).fixed = fp_fe;
        out.(char(A.name)).(char(m)).anova = fp_an;
        out.(char(A.name)).(char(m)).report = fp_md;
        out.(char(A.name)).(char(m)).fig_dir = fp_fig2;
    end
end
end

function y = normalize_high_low(x)
s = lower(strtrim(string(x)));
s(ismember(s,["high","h","1","true","yes","y","高"])) = "high";
s(ismember(s,["low","l","0","false","no","n","低"])) = "low";
y = repmat("", numel(s), 1);
y(s=="high") = "High";
y(s=="low")  = "Low";
end

function w = normalize_wwr(x)
% Accept numeric or strings like WWR45 / 45
s = strtrim(string(x));
% try extract digits
w = repmat("", numel(s), 1);
for i=1:numel(s)
    si = s(i);
    if strlength(si)==0
        continue;
    end
    if ~isnan(str2double(si))
        w(i) = string(str2double(si));
    else
        tok = regexp(char(si), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            w(i) = string(str2double(tok{1}));
        else
            w(i) = si;
        end
    end
end
end

function c = normalize_complexity(x)
% Map 0/1, C0/C1, Low/High to ComplexityLow/ComplexityHigh
s = lower(strtrim(string(x)));

% numeric encodings
s(ismember(s,["0","c0","low","l","complexitylow"])) = "low";
s(ismember(s,["1","c1","high","h","complexityhigh"])) = "high";

c = repmat("", numel(s), 1);
c(s=="low") = "ComplexityLow";
c(s=="high") = "ComplexityHigh";
end

function write_report(fp_md, analysisName, metric, tag, fp_ready, fp_fe, fp_an, lme, AN, fp_fig)
lines = {};
lines{end+1} = sprintf('# TrialIndex LMM (%s) – %s [%s]', string(analysisName), string(metric), string(tag));
lines{end+1} = '';
lines{end+1} = 'Model:';
lines{end+1} = '```';
lines{end+1} = 'EEG ~ Group*WWR*Complexity + TrialIndex + Group:TrialIndex + (1|Subject)';
lines{end+1} = '```';
lines{end+1} = '';
lines{end+1} = 'Key terms to interpret:';
lines{end+1} = '- **TrialIndex**: overall adaptation trend across trials (1–12)';
lines{end+1} = '- **Group:TrialIndex**: whether adaptation slope differs by group';
lines{end+1} = '';
lines{end+1} = 'Outputs:';
lines{end+1} = sprintf('- analysis-ready long table: `%s`', fp_ready);
if strlength(string(fp_fe))>0; lines{end+1} = sprintf('- fixed effects: `%s`', fp_fe); end
if strlength(string(fp_an))>0; lines{end+1} = sprintf('- ANOVA: `%s`', fp_an); end
try
    if nargin>=10 && strlength(string(fp_fig))>0
        lines{end+1} = sprintf('- figures dir: `%s`', fp_fig);
    end
catch
end
lines{end+1} = '';

% Try extract specific coefficient rows
try
    C = lme.Coefficients;
    rn = string(C.Name);
    idxTI = rn=="TrialIndex";
    if any(idxTI)
        lines{end+1} = sprintf('TrialIndex coefficient: Estimate=%.4g, SE=%.4g, t=%.3f, p=%.3g', ...
            C.Estimate(find(idxTI,1)), C.SE(find(idxTI,1)), C.tStat(find(idxTI,1)), C.pValue(find(idxTI,1)));
    end
    % Group:TrialIndex may appear as Group_High:TrialIndex depending on coding
    idxGTI = contains(rn, "TrialIndex") & contains(rn, "Group");
    if any(idxGTI)
        j = find(idxGTI, 1);
        lines{end+1} = sprintf('Group×TrialIndex coefficient (%s): Estimate=%.4g, SE=%.4g, t=%.3f, p=%.3g', ...
            rn(j), C.Estimate(j), C.SE(j), C.tStat(j), C.pValue(j));
    end
catch
end

% ANOVA summary if available
try
    if ~isempty(AN) && height(AN)>0
        % find TrialIndex term row (may be named TrialIndex)
        tt = string(AN.Term);
        i = find(tt=="TrialIndex",1);
        if ~isempty(i)
            lines{end+1} = sprintf('ANOVA TrialIndex: F=%.3f, p=%.3g', AN.FStat(i), AN.pValue(i));
        end
        j = find(contains(tt,"Group") & contains(tt,"TrialIndex"),1);
        if ~isempty(j)
            lines{end+1} = sprintf('ANOVA Group×TrialIndex: F=%.3f, p=%.3g', AN.FStat(j), AN.pValue(j));
        end
    end
catch
end

fid = fopen(fp_md,'w');
if fid<0
    warning('analyze_trialindex_lmm: cannot write %s', fp_md);
    return;
end
fwrite(fid, strjoin(string(lines), newline));
fclose(fid);
end
