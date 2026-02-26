function out = analyze_core_lmm_suite(AllScene, fp_sum, cfg, tag)
%ANALYZE_CORE_LMM_SUITE Task4: Core LMM suite (main effects, 2-way, selective 3-way) + WWR trend checks.
%
% Outputs under:
%   <summary>/analysis-2/task4_core_lmm_suite/
%     factor_WWR/...
%     trend_WWR/...
%
% For each DV (default: O_theta,F_theta,O_alpha,O_beta) and each group type:
%   - Experience (High/Low)
%   - SportFreq  (High/Low)
%
% factor_WWR suite:
%   M1: EEG ~ WWR + Complexity + Group + (1|Subject)
%   M2: EEG ~ WWR*Complexity + WWR*Group + Complexity*Group + (1|Subject)
%   M3: EEG ~ WWR*Complexity*Group + (1|Subject)  (default only for O_alpha,F_theta)
%
% trend_WWR suite (screening / robustness):
%   L:  EEG ~ WWRc + Complexity + Group + (1|Subject)
%   Q:  EEG ~ WWRc + WWRc2 + Complexity + Group + (1|Subject)
%
% Notes:
% - Requires fitlme (Statistics and Machine Learning Toolbox)
% - Requires design columns: WWR, Complexity, Experience/SportFreq

if nargin < 4 || isempty(tag)
    tag = 'raw';
end
if nargin < 3
    cfg = struct();
end

out = struct();

% Check fitlme
if exist('fitlme','file') ~= 2
    warning('analyze_core_lmm_suite: fitlme not found (Stats toolbox missing). Skipping.');
    return;
end

% Required base columns
reqBase = {'subject_id','WWR','Complexity'};
for i=1:numel(reqBase)
    if ~ismember(reqBase{i}, AllScene.Properties.VariableNames)
        warning('analyze_core_lmm_suite: missing required column %s. Skipping.', reqBase{i});
        return;
    end
end

% Metrics
metrics = {"O_theta","F_theta","O_alpha","O_beta"};
try
    if isfield(cfg,'task4_metrics') && ~isempty(cfg.task4_metrics)
        metrics = cellstr(string(cfg.task4_metrics));
    end
catch
end

% Which metrics get 3-way (factor suite)
metrics3 = ["O_alpha","F_theta"];
try
    if isfield(cfg,'task4_factor_threeway_metrics') && ~isempty(cfg.task4_factor_threeway_metrics)
        metrics3 = string(cfg.task4_factor_threeway_metrics);
    end
catch
end

% Output roots
fp_root = fullfile(fp_sum, 'analysis-2', 'task4_core_lmm_suite');
fp_factor = fullfile(fp_root, 'factor_WWR');
fp_trend  = fullfile(fp_root, 'trend_WWR');

% Grouping analyses
analyses = { ...
    struct('name','experience','gcol','ExperienceGroup'), ...
    struct('name','sportfreq','gcol','SportFreqGroup') ...
};

T0 = AllScene;

% normalize columns
T0.subject_id = string(T0.subject_id);
T0.WWR = normalize_wwr(T0.WWR);
T0.Complexity = normalize_complexity(T0.Complexity);

% helper to build dirs
mk = @(p) (exist(p,'dir')||mkdir(p));

for mi = 1:numel(metrics)
    dv = string(metrics{mi});
    if ~ismember(dv, T0.Properties.VariableNames)
        continue;
    end

    for ai=1:numel(analyses)
        A = analyses{ai};
        gcol = string(A.gcol);
        if ~ismember(gcol, string(T0.Properties.VariableNames))
            continue;
        end
        gcolc = char(gcol);
        dvc = char(dv);

        % Prepare analysis table
        T = T0(:, {'subject_id','WWR','Complexity',gcolc,dvc});
        T.Properties.VariableNames{4} = 'GroupRaw';
        % GroupRaw comes from canonical group columns (ExperienceGroup/SportFreqGroup)
        T.EEG = double(T.(dvc));
        T.Group = normalize_high_low(T.GroupRaw);
        T = remove_missing_rows(T);
        if height(T) < 20
            continue;
        end

        % Ensure categorical factors
        T.Subject = categorical(string(T.subject_id));
        T.WWR = categorical(string(T.WWR), {'15','45','75'});
        T.Complexity = categorical(string(T.Complexity), {'ComplexityLow','ComplexityHigh'});
        T.Group = categorical(string(T.Group), {'Low','High'});

        % ===== factor suite =====
        fp_tbl = fullfile(fp_factor, 'tables', tag, A.name);
        fp_rep = fullfile(fp_factor, 'reports', tag, A.name);
        fp_fig = fullfile(fp_factor, 'figures', tag, A.name);
        mk(fp_tbl); mk(fp_rep); mk(fp_fig);

        % Save analysis-ready
        fp_ready = fullfile(fp_tbl, sprintf('analysis_ready_%s_%s.csv', dv, tag));
        try
            writetable(T(:,{'Subject','WWR','Complexity','Group','EEG'}), fp_ready);
        catch
        end

        % Fit models
        try
            lme1 = fitlme(T, 'EEG ~ WWR + Complexity + Group + (1|Subject)');
            lme2 = fitlme(T, 'EEG ~ WWR*Complexity + WWR*Group + Complexity*Group + (1|Subject)');
        catch ME
            warning('analyze_core_lmm_suite: fitlme failed (%s/%s/%s): %s', tag, A.name, dv, ME.message);
            continue;
        end

        write_lme_tables(lme1, fp_tbl, sprintf('model1_main_effects_%s_%s', dv, tag));
        write_lme_tables(lme2, fp_tbl, sprintf('model2_two_way_%s_%s', dv, tag));

        % selective 3-way
        did3 = false;
        lme3 = [];
        if any(dv == metrics3)
            try
                lme3 = fitlme(T, 'EEG ~ WWR*Complexity*Group + (1|Subject)');
                write_lme_tables(lme3, fp_tbl, sprintf('model3_three_way_%s_%s', dv, tag));
                did3 = true;
            catch ME
                warning('analyze_core_lmm_suite: model3 failed (%s/%s/%s): %s', tag, A.name, dv, ME.message);
            end
        end

        % Simple direction summaries (raw means)
        Summ = summarize_direction(T);
        try
            writetable(Summ, fullfile(fp_tbl, sprintf('direction_means_%s_%s.csv', dv, tag)));
        catch
        end

        % Figures (paper-friendly basics)
        try
            plot_factor_figures(T, dv, tag, A.name, fp_fig);
        catch ME
            fprintf(2,'[WARN] task4 plot_factor_figures failed for %s/%s/%s: %s\n', tag, A.name, dv, ME.message);
        end

        % Report
        try
            fp_md = fullfile(fp_rep, sprintf('task4_factor_report_%s_%s.md', dv, tag));
            write_factor_report(fp_md, dv, tag, A, lme1, lme2, lme3, did3, Summ);
        catch
        end

        % ===== trend suite =====
        fp_tbl_t = fullfile(fp_trend, 'tables', tag, A.name);
        fp_rep_t = fullfile(fp_trend, 'reports', tag, A.name);
        fp_fig_t = fullfile(fp_trend, 'figures', tag, A.name);
        mk(fp_tbl_t); mk(fp_rep_t); mk(fp_fig_t);

        % Build numeric WWRc
        w = double(string(T.WWR));
        w = str2double(string(T.WWR)); %#ok<NASGU>
        % T.WWR is categorical now; recover from categories
        wnum = str2double(string(T.WWR));
        WWRc = wnum - 45; % center at 45
        WWRc2 = WWRc.^2;
        Tt = table();
        Tt.Subject = T.Subject;
        Tt.Complexity = T.Complexity;
        Tt.Group = T.Group;
        Tt.WWRc = WWRc;
        Tt.WWRc2 = WWRc2;
        Tt.EEG = T.EEG;

        try
            lmeL = fitlme(Tt, 'EEG ~ WWRc + Complexity + Group + (1|Subject)');
            lmeQ = fitlme(Tt, 'EEG ~ WWRc + WWRc2 + Complexity + Group + (1|Subject)');
        catch ME
            warning('analyze_core_lmm_suite: trend fitlme failed (%s/%s/%s): %s', tag, A.name, dv, ME.message);
            continue;
        end

        write_lme_tables(lmeL, fp_tbl_t, sprintf('trend_linear_%s_%s', dv, tag));
        write_lme_tables(lmeQ, fp_tbl_t, sprintf('trend_quadratic_%s_%s', dv, tag));

        try
            fp_md_t = fullfile(fp_rep_t, sprintf('task4_trend_report_%s_%s.md', dv, tag));
            write_trend_report(fp_md_t, dv, tag, A, lmeL, lmeQ);
        catch
        end

        try
            plot_trend_figures(Tt, dv, tag, A.name, fp_fig_t);
        catch
        end

    end
end

end

% ===== helpers =====

function T = remove_missing_rows(T)
T = T(~isnan(T.EEG),:);
T = T(strlength(string(T.WWR))>0,:);
T = T(strlength(string(T.Complexity))>0,:);
T = T(strlength(string(T.Group))>0,:);
T = T(string(T.Group)=="High" | string(T.Group)=="Low",:);
end

function write_lme_tables(lme, fp_tbl, stem)
% fixed effects
try
    C = lme.Coefficients;
    writetable(C, fullfile(fp_tbl, sprintf('%s_fixed_effects.csv', stem)));
catch
end
% anova table
try
    A = anova(lme,'DFMethod','Satterthwaite');
    writetable(A, fullfile(fp_tbl, sprintf('%s_anova.csv', stem)));
catch
end
end

function Summ = summarize_direction(T)
% Raw mean±SEM by factor levels (not model-adjusted)
Summ = table();
try
    % WWR
    [g, w] = findgroups(T.WWR);
    mu = splitapply(@(x) mean(x,'omitnan'), T.EEG, g);
    se = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T.EEG, g);
    n = splitapply(@numel, T.EEG, g);
    Tw = table(string(w), n, mu, se, 'VariableNames',{'level','n','mean','sem'});
    Tw.factor = repmat("WWR",height(Tw),1);

    % Complexity
    [g, c] = findgroups(T.Complexity);
    mu = splitapply(@(x) mean(x,'omitnan'), T.EEG, g);
    se = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T.EEG, g);
    n = splitapply(@numel, T.EEG, g);
    Tc = table(string(c), n, mu, se, 'VariableNames',{'level','n','mean','sem'});
    Tc.factor = repmat("Complexity",height(Tc),1);

    % Group
    [g, gr] = findgroups(T.Group);
    mu = splitapply(@(x) mean(x,'omitnan'), T.EEG, g);
    se = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T.EEG, g);
    n = splitapply(@numel, T.EEG, g);
    Tg = table(string(gr), n, mu, se, 'VariableNames',{'level','n','mean','sem'});
    Tg.factor = repmat("Group",height(Tg),1);

    Summ = [Tw; Tc; Tg];
catch
end
end

function plot_factor_figures(T, dv, tag, analysisName, fp_fig)
set(0,'DefaultFigureVisible','off');

% 1) mean by WWR
fig = figure('Color','w','Position',[60 60 850 320]);
ax = subplot(1,3,1); hold(ax,'on');
plot_mean_by(ax, T, 'WWR'); title(ax,'Mean by WWR');
ax2 = subplot(1,3,2); hold(ax2,'on');
plot_lines(ax2, T, 'WWR', 'Complexity'); title(ax2,'WWR × Complexity');
ax3 = subplot(1,3,3); hold(ax3,'on');
plot_lines(ax3, T, 'WWR', 'Group'); title(ax3,'WWR × Group');
sgtitle(sprintf('Task4 factor | %s | %s [%s]', analysisName, dv, tag), 'Interpreter','none');
fp = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('task4_factor_%s_%s_%s.png', tag, analysisName, dv)));
pipeline.export_figure_png(fig, fp, 300);
try; close(fig); catch; end
end

function plot_mean_by(ax, T, fac)
trials = categories(T.(fac));
mu = zeros(numel(trials),1); se=mu;
for i=1:numel(trials)
    yy = T.EEG(T.(fac)==trials{i});
    mu(i)=mean(yy,'omitnan');
    se(i)=std(yy,'omitnan')/sqrt(numel(yy));
end
errorbar(ax, 1:numel(trials), mu, se, 'o-','LineWidth',1.4);
set(ax,'XTick',1:numel(trials),'XTickLabel',trials);
xtickangle(ax,0);
grid(ax,'on');
xlabel(ax, fac);
end

function plot_lines(ax, T, xfac, linefac)
xl = categories(T.(xfac));
ll = categories(T.(linefac));
cols = lines(numel(ll));
for j=1:numel(ll)
    mu = nan(numel(xl),1); se=mu;
    for i=1:numel(xl)
        yy = T.EEG(T.(xfac)==xl{i} & T.(linefac)==ll{j});
        if isempty(yy); continue; end
        mu(i)=mean(yy,'omitnan');
        se(i)=std(yy,'omitnan')/sqrt(numel(yy));
    end
    errorbar(ax, 1:numel(xl), mu, se, 'o-','LineWidth',1.2,'Color',cols(j,:),...
        'MarkerFaceColor',cols(j,:), 'DisplayName', ll{j});
end
set(ax,'XTick',1:numel(xl),'XTickLabel',xl);
grid(ax,'on');
xlabel(ax, xfac);
legend(ax,'Location','best');
end

function write_factor_report(fp_md, dv, tag, A, lme1, lme2, lme3, did3, Summ)
lines = strings(0,1);
lines(end+1) = sprintf('# Task4 Core LMM (factor WWR) | %s | %s | %s', A.name, dv, tag);
lines(end+1) = '';
lines(end+1) = '## Models';
lines(end+1) = '- Model1: `EEG ~ WWR + Complexity + Group + (1|Subject)`';
lines(end+1) = '- Model2: `EEG ~ WWR*Complexity + WWR*Group + Complexity*Group + (1|Subject)`';
if did3
    lines(end+1) = '- Model3: `EEG ~ WWR*Complexity*Group + (1|Subject)` (selective)';
end
lines(end+1) = '';
lines(end+1) = '## Quick direction (raw means, not model-adjusted)';
try
    % summarize as markdown table
    lines = [lines; ""; "factor | level | n | mean | sem"; "---|---:|---:|---:|---:"];
    for i=1:height(Summ)
        lines(end+1) = sprintf('%s|%s|%d|%.6g|%.6g', string(Summ.factor(i)), string(Summ.level(i)), Summ.n(i), Summ.mean(i), Summ.sem(i));
    end
catch
end
lines(end+1) = '';
lines(end+1) = '## ANOVA (Satterthwaite)';
lines(end+1) = 'See exported CSVs in tables/.';
lines(end+1) = '';
lines(end+1) = '## Notes';
lines(end+1) = '- Reference levels: WWR=15, Complexity=Low, Group=Low (categorical order enforced).';
lines(end+1) = '- Use Model1 for main effects, Model2 for two-way interactions; only interpret three-way when fitted.';

fid = fopen(fp_md,'w');
for i=1:numel(lines)
    fprintf(fid,'%s\n', lines(i));
end
fclose(fid);
end

function plot_trend_figures(Tt, dv, tag, analysisName, fp_fig)
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[80 80 740 320]);
ax = axes(fig); hold(ax,'on');
% plot raw means by WWR value
wwrVals = [-30 0 30];
labels = ["15","45","75"];
mu = nan(3,1); se=mu;
for i=1:3
    yy = Tt.EEG(Tt.WWRc==wwrVals(i));
    mu(i)=mean(yy,'omitnan');
    se(i)=std(yy,'omitnan')/sqrt(numel(yy));
end
errorbar(ax, 1:3, mu, se, 'o-','LineWidth',1.4);
set(ax,'XTick',1:3,'XTickLabel',labels);
xlabel(ax,'WWR'); ylabel(ax, dv);
grid(ax,'on');
title(ax, sprintf('Task4 trend screen | %s | %s [%s]', analysisName, dv, tag),'Interpreter','none');
fp = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('task4_trend_%s_%s_%s.png', tag, analysisName, dv)));
pipeline.export_figure_png(fig, fp, 300);
try; close(fig); catch; end
end

function write_trend_report(fp_md_t, dv, tag, A, lmeL, lmeQ)
lines = strings(0,1);
lines(end+1) = sprintf('# Task4 WWR Trend Screen | %s | %s | %s', A.name, dv, tag);
lines(end+1) = '';
lines(end+1) = '## Models';
lines(end+1) = '- Linear: `EEG ~ WWRc + Complexity + Group + (1|Subject)`';
lines(end+1) = '- Quadratic: `EEG ~ WWRc + WWRc2 + Complexity + Group + (1|Subject)`';
lines(end+1) = '';
lines(end+1) = '## How to read';
lines(end+1) = '- `WWRc` tests monotonic trend (increasing/decreasing with WWR).';
lines(end+1) = '- `WWRc2` tests curvature (U/∩ shape). With only 3 levels, treat as screening evidence; confirm with PeakIndex task.';
lines(end+1) = '';
lines(end+1) = '## Outputs';
lines(end+1) = 'See exported fixed_effects/anova CSVs in tables/.';

fid = fopen(fp_md_t,'w');
for i=1:numel(lines)
    fprintf(fid,'%s\n', lines(i));
end
fclose(fid);
end

function w = normalize_wwr(x)
% normalize to strings {15,45,75}
w = string(x);
w = strtrim(w);
% extract digits if needed
for i=1:numel(w)
    if strlength(w(i))==0; continue; end
    tok = regexp(char(w(i)), '(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        w(i) = string(str2double(tok{1}));
    end
end
% keep only supported
ok = ismember(w,["15","45","75"]);
w(~ok) = "";
end

function c = normalize_complexity(x)
% normalize to ComplexityLow/High
c = string(x);
c = strtrim(c);
cl = lower(c);
out = repmat("", numel(c), 1);
out(ismember(cl,["low","0","c0","complexitylow"])) = "ComplexityLow";
out(ismember(cl,["high","1","c1","complexityhigh"])) = "ComplexityHigh";
% if already canonical
out(c=="ComplexityLow") = "ComplexityLow";
out(c=="ComplexityHigh") = "ComplexityHigh";
% sometimes numeric
isNum = ~isnan(str2double(c));
out(isNum & str2double(c)==0) = "ComplexityLow";
out(isNum & str2double(c)==1) = "ComplexityHigh";
c = out;
end

function g = normalize_high_low(x)
% normalize to High/Low
s = string(x);
s = strtrim(s);
sl = lower(s);
g = repmat("", numel(s), 1);
g(ismember(sl,["high","1","高","h"])) = "High";
g(ismember(sl,["low","0","低","l"])) = "Low";
% already
mask = (s=="High" | s=="Low");
g(mask) = s(mask);
end
