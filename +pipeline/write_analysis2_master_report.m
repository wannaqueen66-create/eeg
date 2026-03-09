function out = write_analysis2_master_report(fp_sum, cfg)
%WRITE_ANALYSIS2_MASTER_REPORT Build one ordered master report for analysis-2.
%
% Output:
%   <summary>/analysis-2/MASTER_REPORT.md
%
% Report order (requested):
% 1) Main effects
% 2) Two-way interactions
% 3) Three-way interactions (only when significant)
% 4) Sequence control (TrialIndex)

if nargin < 2
    cfg = struct();
end

out = struct();
if exist(fullfile(fp_sum, 'analysis'),'dir')
    fp_a2 = fullfile(fp_sum, 'analysis');
else
    fp_a2 = fullfile(fp_sum, 'analysis-2');
end
if ~exist(fp_a2,'dir')
    return;
end

tags = {"raw","qc"};
try
    if isfield(cfg,'analysis2_master_tags') && ~isempty(cfg.analysis2_master_tags)
        tags = cellstr(string(cfg.analysis2_master_tags));
    end
catch
end

analyses = {"experience","sportfreq"};
metrics = {"O_theta","F_theta","O_alpha","O_beta"};

lines = {};
add('##HEADER##');
add('# Analysis MASTER REPORT');
add('');
add(sprintf('Generated: %s', datestr(now,31)));
add('');
add('This file enforces the requested reporting order:');
add('1) Main effects  2) Two-way interactions  3) Three-way (only significant)  4) Sequence control');
add('');

for it = 1:numel(tags)
    tag = string(tags{it});
    add('---');
    add(sprintf('## Data tag: `%s`', tag));
    add('');

    for ia = 1:numel(analyses)
        a = string(analyses{ia});
        add(sprintf('### Grouping: `%s`', a));
        add('');

        % 1) main effects
        add('#### 1) Main effects (Task4 Model1)');
        rows = collect_task4_terms(fp_sum, tag, a, metrics, 'model1_main_effects', ["WWR","Complexity","Group"]);
        append_lines(render_rows(rows, true));
        add('');

        % 2) two-way
        add('#### 2) Two-way interactions (Task4 Model2)');
        rows2 = collect_task4_terms(fp_sum, tag, a, metrics, 'model2_two_way', ["WWR:Complexity","WWR:Group","Complexity:Group"]);
        append_lines(render_rows(rows2, false));
        add('');

        % 3) three-way only significant
        add('#### 3) Three-way interaction (Task4 Model3, only significant shown)');
        rows3 = collect_task4_terms(fp_sum, tag, a, metrics, 'model3_three_way', ["WWR:Complexity:Group"]);
        if isempty(rows3)
            add('- no model3 files found');
        else
            p = toNumField(rows3,'p');
            keep = ~isnan(p) & p < 0.05;
            if any(keep)
                append_lines(render_rows(rows3(keep,:), false));
            else
                add('- no significant three-way interaction (p < 0.05)');
            end
        end
        add('');

        % 4) sequence control
        add('#### 4) Sequence control (Task3 TrialIndex model)');
        rows4 = collect_task3_terms(fp_sum, tag, a, metrics, ["TrialIndex","Group:TrialIndex"]);
        append_lines(render_rows(rows4, false));
        add('');

        % Significant-only digest + one-line direction summary
        add('#### Significant-only digest (p < 0.05)');
        Sig = build_sig_digest(rows, rows2, rows3, rows4);
        append_lines(render_sig_rows(Sig));
        add(summarize_sig_sentence(Sig, tag, a));
        add('');

        % links
        add('- Task5 PeakIndex files: `analysis/task5_peakindex_invertedu/...` (legacy fallback may still appear under `analysis-2/...`)');
        add('- Task6 core-metric special files: `analysis/task6_coremetric_special/...` (legacy fallback may still appear under `analysis-2/task6_obeta_special/...`)');
        add('- Task7 individual checks: `analysis/task7_individual_checks/...` (legacy fallback may still appear under `analysis-2/...`)');
        add('');
    end
end

fp = fullfile(fp_a2, 'MASTER_REPORT.md');
fid = fopen(fp,'w');
for i=1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
fclose(fid);
out.master_report = fp;

    function add(s)
        if strcmp(s,'##HEADER##'); return; end
        lines{end+1} = char(string(s));
    end

    function append_lines(c)
        if isempty(c); return; end
        for ii=1:numel(c)
            lines{end+1} = char(string(c{ii}));
        end
    end
end

function rows = collect_task4_terms(fp_sum, tag, analysisName, metrics, stemPrefix, wantedTerms)
rows = table();
if exist(fullfile(fp_sum,'analysis'),'dir')
    base = fullfile(fp_sum,'analysis','task4_core_lmm_suite',char(tag),'tables','factor_WWR',char(analysisName));
else
    base = fullfile(fp_sum,'analysis-2','task4_core_lmm_suite','factor_WWR','tables',char(tag),char(analysisName));
end
if ~exist(base,'dir'); return; end

for im = 1:numel(metrics)
    m = string(metrics{im});
    fpAn = fullfile(base, sprintf('%s_%s_%s_anova.csv', stemPrefix, m, tag));
    if ~exist(fpAn,'file'); continue; end
    T = safe_read(fpAn);
    if isempty(T); continue; end

    [termCol,pCol,fCol] = find_cols(T);
    if termCol=="" || pCol==""; continue; end
    terms = string(T.(termCol));
    p = toDouble(T.(pCol));
    F = nan(size(p));
    if fCol ~= ""
        F = toDouble(T.(fCol));
    end

    for wi=1:numel(wantedTerms)
        wt = wantedTerms(wi);
        idx = find(contains(terms, wt), 1, 'first');
        if isempty(idx); continue; end
        direction = "";
        if wt=="Group"
            direction = read_group_direction(fp_sum, tag, analysisName, m);
        end
        rows = [rows; table(m, wt, F(idx), p(idx), direction, ...
            'VariableNames', {'metric','term','F','p','direction'})]; %#ok<AGROW>
    end
end
end

function rows = collect_task3_terms(fp_sum, tag, analysisName, metrics, wantedTerms)
rows = table();
if exist(fullfile(fp_sum,'analysis'),'dir')
    base = fullfile(fp_sum,'analysis','task3_trialindex_lmm',char(tag),'tables',char(analysisName));
else
    base = fullfile(fp_sum,'analysis-2','task3_trialindex_lmm','tables',char(tag),char(analysisName));
end
if ~exist(base,'dir'); return; end

for im = 1:numel(metrics)
    m = string(metrics{im});
    fpAn = fullfile(base, sprintf('lmm_anova_%s_%s.csv', m, tag));
    if ~exist(fpAn,'file'); continue; end
    T = safe_read(fpAn);
    if isempty(T); continue; end

    [termCol,pCol,fCol] = find_cols(T);
    if termCol=="" || pCol==""; continue; end
    terms = string(T.(termCol));
    p = toDouble(T.(pCol));
    F = nan(size(p));
    if fCol ~= ""
        F = toDouble(T.(fCol));
    end

    for wi=1:numel(wantedTerms)
        wt = wantedTerms(wi);
        idx = find(contains(terms, wt), 1, 'first');
        if isempty(idx) && wt=="Group:TrialIndex"
            idx = find(contains(terms,'TrialIndex:Group'),1,'first');
        end
        if isempty(idx); continue; end
        rows = [rows; table(m, wt, F(idx), p(idx), ...
            'VariableNames', {'metric','term','F','p'})]; %#ok<AGROW>
    end
end
end

function dirText = read_group_direction(fp_sum, tag, analysisName, metric)
dirText = "";
if exist(fullfile(fp_sum,'analysis'),'dir')
    fp = fullfile(fp_sum,'analysis','task4_core_lmm_suite',char(tag),'tables','factor_WWR',char(analysisName), ...
        sprintf('direction_means_%s_%s.csv', metric, tag));
else
    fp = fullfile(fp_sum,'analysis-2','task4_core_lmm_suite','factor_WWR','tables',char(tag),char(analysisName), ...
        sprintf('direction_means_%s_%s.csv', metric, tag));
end
if ~exist(fp,'file'); return; end
T = safe_read(fp);
if isempty(T); return; end
if ~all(ismember({'factor','level','mean'}, string(T.Properties.VariableNames)))
    return;
end
G = T(string(T.factor)=="Group", :);
if height(G) < 2; return; end
mh = NaN; ml = NaN;
try
    mh = double(G.mean(string(G.level)=="High"));
    ml = double(G.mean(string(G.level)=="Low"));
catch
end
if isempty(mh) || isempty(ml); return; end
if mh > ml
    dirText = "High > Low";
elseif mh < ml
    dirText = "Low > High";
else
    dirText = "High ~= Low";
end
end

function T = safe_read(fp)
try
    T = readtable(fp, 'TextType','string');
catch
    T = table();
end
end

function [termCol,pCol,fCol] = find_cols(T)
vars = string(T.Properties.VariableNames);
termCol = first_match(vars, ["Term","Name","Effect","Source"]);
pCol = first_match(vars, ["pValue","pvalue","p","ProbF","PrF"]);
fCol = first_match(vars, ["FStat","F","FValue","F_stat","Fstat"]);
end

function c = first_match(vars, cands)
c = "";
for i=1:numel(cands)
    if any(vars==cands(i))
        c = cands(i);
        return;
    end
end
end

function x = toDouble(v)
try
    x = double(v);
catch
    try
        x = str2double(string(v));
    catch
        x = nan(size(v));
    end
end
end

function x = toNumField(T, name)
if ismember(name, string(T.Properties.VariableNames))
    x = toDouble(T.(name));
else
    x = nan(height(T),1);
end
end

function c = render_rows(T, withDirection)
c = {};
if isempty(T) || height(T)==0
    c{end+1} = '- no rows found';
    return;
end
if withDirection && ismember('direction', string(T.Properties.VariableNames))
    c{end+1} = 'metric | term | F | p | direction';
    c{end+1} = '---|---|---:|---:|---';
    for i=1:height(T)
        c{end+1} = sprintf('%s | %s | %.4g | %.4g | %s', string(T.metric(i)), string(T.term(i)), double(T.F(i)), double(T.p(i)), string(T.direction(i)));
    end
else
    c{end+1} = 'metric | term | F | p';
    c{end+1} = '---|---|---:|---:';
    for i=1:height(T)
        c{end+1} = sprintf('%s | %s | %.4g | %.4g', string(T.metric(i)), string(T.term(i)), double(T.F(i)), double(T.p(i)));
    end
end
end

function Sig = build_sig_digest(rowsMain, rows2way, rows3way, rowsSeq)
Sig = table('Size',[0 7], ...
    'VariableTypes', {'string','string','string','double','double','string','string'}, ...
    'VariableNames', {'section','metric','term','F','p','direction','note'});

if ~isempty(rowsMain)
    p = toNumField(rowsMain,'p');
    k = find(~isnan(p) & p<0.05);
    for i=k'
        note = "";
        if rowsMain.term(i)=="Group"
            note = "main-effect group direction";
        end
        Sig = [Sig; {"main", string(rowsMain.metric(i)), string(rowsMain.term(i)), double(rowsMain.F(i)), double(rowsMain.p(i)), safe_dir(rowsMain,i), note}]; %#ok<AGROW>
    end
end
if ~isempty(rows2way)
    p = toNumField(rows2way,'p');
    k = find(~isnan(p) & p<0.05);
    for i=k'
        Sig = [Sig; {"two-way", string(rows2way.metric(i)), string(rows2way.term(i)), double(rows2way.F(i)), double(rows2way.p(i)), "", "interaction significant"}]; %#ok<AGROW>
    end
end
if ~isempty(rows3way)
    p = toNumField(rows3way,'p');
    k = find(~isnan(p) & p<0.05);
    for i=k'
        Sig = [Sig; {"three-way", string(rows3way.metric(i)), string(rows3way.term(i)), double(rows3way.F(i)), double(rows3way.p(i)), "", "key higher-order interaction"}]; %#ok<AGROW>
    end
end
if ~isempty(rowsSeq)
    p = toNumField(rowsSeq,'p');
    k = find(~isnan(p) & p<0.05);
    for i=k'
        nt = "";
        if rowsSeq.term(i)=="TrialIndex"
            nt = "adaptation trend present";
        elseif rowsSeq.term(i)=="Group:TrialIndex"
            nt = "different adaptation slopes by group";
        end
        Sig = [Sig; {"sequence", string(rowsSeq.metric(i)), string(rowsSeq.term(i)), double(rowsSeq.F(i)), double(rowsSeq.p(i)), "", nt}]; %#ok<AGROW>
    end
end
end

function d = safe_dir(T, i)
d = "";
try
    if ismember('direction', string(T.Properties.VariableNames))
        d = string(T.direction(i));
    end
catch
end
end

function c = render_sig_rows(Sig)
c = {};
if isempty(Sig) || height(Sig)==0
    c{end+1} = '- no significant effects (p < 0.05) for this grouping/tag';
    return;
end
c{end+1} = 'section | metric | term | F | p | direction | note';
c{end+1} = '---|---|---|---:|---:|---|---';
for i=1:height(Sig)
    c{end+1} = sprintf('%s | %s | %s | %.4g | %.4g | %s | %s', ...
        string(Sig.section(i)), string(Sig.metric(i)), string(Sig.term(i)), ...
        double(Sig.F(i)), double(Sig.p(i)), string(Sig.direction(i)), string(Sig.note(i)));
end
end

function s = summarize_sig_sentence(Sig, tag, analysisName)
if isempty(Sig) || height(Sig)==0
    s = sprintf('- Summary (%s/%s): no significant effects detected at p < 0.05 in the extracted core terms.', tag, analysisName);
    return;
end
nMain = sum(string(Sig.section)=="main");
n2 = sum(string(Sig.section)=="two-way");
n3 = sum(string(Sig.section)=="three-way");
nS = sum(string(Sig.section)=="sequence");
idxDir = find(strlength(string(Sig.direction))>0);
dirTxt = "";
if ~isempty(idxDir)
    d = unique(string(Sig.direction(idxDir)), 'stable');
    dirTxt = strjoin(d, '; ');
end
if strlength(dirTxt)>0
    s = sprintf('- Summary (%s/%s): sig main=%d, two-way=%d, three-way=%d, sequence=%d. Group direction notes: %s.', ...
        tag, analysisName, nMain, n2, n3, nS, dirTxt);
else
    s = sprintf('- Summary (%s/%s): sig main=%d, two-way=%d, three-way=%d, sequence=%d.', ...
        tag, analysisName, nMain, n2, n3, nS);
end
end
