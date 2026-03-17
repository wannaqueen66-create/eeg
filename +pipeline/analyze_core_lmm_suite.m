function out = analyze_core_lmm_suite(AllScene, fp_sum, cfg, tag)
%ANALYZE_CORE_LMM_SUITE Task4: Core LMM suite (main effects, 2-way, selective 3-way) + WWR trend checks.
% Also optionally runs factor-suite split by round/block (Block1 vs Block2) when block_id exists.
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
% Default: run 3-way on all core metrics (can still override via config)
metrics3 = ["O_theta","F_theta","O_alpha","O_beta"];
try
    if isfield(cfg,'task4_factor_threeway_metrics') && ~isempty(cfg.task4_factor_threeway_metrics)
        metrics3 = string(cfg.task4_factor_threeway_metrics);
    end
catch
end

% Output roots
[fp_root_base, fp_task_tbl, fp_task_fig, fp_task_rep] = pipeline.get_analysis_task_subdirs(fp_sum, 'task4_core_lmm_suite', tag);
fp_root = fp_root_base;
fp_factor = fullfile(fp_root, 'factor_WWR');
fp_trend  = fullfile(fp_root, 'trend_WWR');

% Grouping analyses
analyses = { ...
    struct('name','experience','gcol','ExperienceGroup') ...
};

T0 = AllScene;

% normalize columns
T0.subject_id = string(T0.subject_id);
T0.WWR = normalize_wwr(T0.WWR);
T0.Complexity = normalize_complexity(T0.Complexity);

% helper to build dirs
mk = @(p) (exist(p,'dir')||mkdir(p));

FactorSummary = table();
FactorRoundSummary = table();
TrendSummary = table();
TrendRoundSummary = table();

% Optional: also run factor-suite split by round/block when block_id exists.
% Default: true (requested), can disable via cfg.task4_factor_split_by_round = false.
doFactorRoundSplit = true;
try
    if isfield(cfg,'task4_factor_split_by_round')
        doFactorRoundSplit = logical(cfg.task4_factor_split_by_round);
    end
catch
end

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
        vars = {'subject_id','WWR','Complexity',gcolc,dvc};
        if ismember('block_id', T0.Properties.VariableNames)
            vars{end+1} = 'block_id';
        end
        T = T0(:, vars);
        T.Properties.VariableNames{4} = 'GroupRaw';
        % GroupRaw comes from canonical group columns (ExperienceGroup)
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

        % ===== factor suite (overall) =====
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

        % Collect factor summary row (paper overview)
        try
            [pWWR, pCx, pGrp, pWxC, pWxG, pCxG, p3] = extract_factor_pvals(lme1, lme2, lme3, did3);
            rr = table(string(A.name), dv, "overall", NaN, pWWR, pCx, pGrp, pWxC, pWxG, pCxG, p3, ...
                'VariableNames', {'analysis','metric','scope','round','p_WWR','p_Complexity','p_Group','p_WWRxComplexity','p_WWRxGroup','p_ComplexityxGroup','p_threeway'});
            FactorSummary = [FactorSummary; rr]; %#ok<AGROW>
        catch
        end

        % ===== factor suite (split by round/block) =====
        try
            if doFactorRoundSplit && ismember('block_id', T.Properties.VariableNames)
                rounds = unique(double(T.block_id));
                rounds = rounds(isfinite(rounds));
                rounds = rounds(:)';
                rounds = rounds(rounds==1 | rounds==2);
                for rb = rounds
                    Tb = T(double(T.block_id)==rb, :);
                    if height(Tb) < 20
                        continue;
                    end

                    fp_tbl_r = fullfile(fp_factor, 'tables', tag, A.name, sprintf('round%d', rb));
                    fp_rep_r = fullfile(fp_factor, 'reports', tag, A.name, sprintf('round%d', rb));
                    fp_fig_r = fullfile(fp_factor, 'figures', tag, A.name, sprintf('round%d', rb));
                    mk(fp_tbl_r); mk(fp_rep_r); mk(fp_fig_r);

                    fp_ready_r = fullfile(fp_tbl_r, sprintf('analysis_ready_%s_%s_round%d.csv', dv, tag, rb));
                    try
                        writetable(Tb(:,{'Subject','WWR','Complexity','Group','EEG'}), fp_ready_r);
                    catch
                    end

                    try
                        lme1b = fitlme(Tb, 'EEG ~ WWR + Complexity + Group + (1|Subject)');
                        lme2b = fitlme(Tb, 'EEG ~ WWR*Complexity + WWR*Group + Complexity*Group + (1|Subject)');
                    catch ME
                        warning('analyze_core_lmm_suite: round factor fitlme failed (%s/%s/%s/r%d): %s', tag, A.name, dv, rb, ME.message);
                        continue;
                    end

                    write_lme_tables(lme1b, fp_tbl_r, sprintf('model1_main_effects_%s_%s_round%d', dv, tag, rb));
                    write_lme_tables(lme2b, fp_tbl_r, sprintf('model2_two_way_%s_%s_round%d', dv, tag, rb));

                    did3b = false;
                    lme3b = [];
                    if any(dv == metrics3)
                        try
                            lme3b = fitlme(Tb, 'EEG ~ WWR*Complexity*Group + (1|Subject)');
                            write_lme_tables(lme3b, fp_tbl_r, sprintf('model3_three_way_%s_%s_round%d', dv, tag, rb));
                            did3b = true;
                        catch ME
                            warning('analyze_core_lmm_suite: round model3 failed (%s/%s/%s/r%d): %s', tag, A.name, dv, rb, ME.message);
                        end
                    end

                    Summb = summarize_direction(Tb);
                    try
                        writetable(Summb, fullfile(fp_tbl_r, sprintf('direction_means_%s_%s_round%d.csv', dv, tag, rb)));
                    catch
                    end

                    try
                        plot_factor_figures(Tb, dv, tag, sprintf('%s_round%d', A.name, rb), fp_fig_r);
                    catch
                    end

                    try
                        fp_mdr = fullfile(fp_rep_r, sprintf('task4_factor_report_%s_%s_round%d.md', dv, tag, rb));
                        write_factor_report(fp_mdr, dv, tag, struct('name',sprintf('%s_round%d',A.name,rb),'gcol',A.gcol), lme1b, lme2b, lme3b, did3b, Summb);
                    catch
                    end

                    try
                        [pWWR, pCx, pGrp, pWxC, pWxG, pCxG, p3] = extract_factor_pvals(lme1b, lme2b, lme3b, did3b);
                        rr = table(string(A.name), dv, "round", rb, pWWR, pCx, pGrp, pWxC, pWxG, pCxG, p3, ...
                            'VariableNames', {'analysis','metric','scope','round','p_WWR','p_Complexity','p_Group','p_WWRxComplexity','p_WWRxGroup','p_ComplexityxGroup','p_threeway'});
                        FactorRoundSummary = [FactorRoundSummary; rr]; %#ok<AGROW>
                    catch
                    end
                end
            end
        catch ME
            fprintf(2,'[WARN] task4 factor round split failed (%s/%s/%s): %s\n', tag, A.name, dv, ME.message);
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
        if ismember('block_id', T.Properties.VariableNames)
            Tt.block_id = double(T.block_id);
        end

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

        % Collect trend summary row (paper overview)
        try
            [bL, pL, bQ2, pQ2] = extract_trend_stats(lmeL, lmeQ);
            rt = table(string(A.name), dv, bL, pL, bQ2, pQ2, ...
                'VariableNames', {'analysis','metric','beta_linear','p_linear','beta_quadratic','p_quadratic'});
            TrendSummary = [TrendSummary; rt]; %#ok<AGROW>
        catch
        end

        % Round-specific (Block1 vs Block2) quadratic checks for inverted-U consistency
        try
            if ismember('block_id', Tt.Properties.VariableNames)
                TR = collect_round_trend_stats(Tt, A.name, dv);
                if ~isempty(TR) && height(TR)>0
                    TrendRoundSummary = [TrendRoundSummary; TR]; %#ok<AGROW>
                    fp_round = fullfile(fp_tbl_t, sprintf('trend_round_summary_%s_%s.csv', dv, tag));
                    writetable(TR, fp_round);
                end
            end
        catch ME
            fprintf(2,'[WARN] task4 round trend failed (%s/%s/%s): %s\n', tag, A.name, dv, ME.message);
        end

    end
end

% Write overview summary tables + paper-friendly overview figures
try
    if ~isempty(FactorSummary) && height(FactorSummary)>0
        fp_tbl_over = fullfile(fp_factor, 'tables', tag, 'overview');
        fp_fig_over = fullfile(fp_factor, 'figures', tag, 'overview');
        mk(fp_tbl_over); mk(fp_fig_over);
        fp_csv = fullfile(fp_tbl_over, sprintf('task4_factor_summary_%s.csv', tag));
        writetable(FactorSummary, fp_csv);
        plot_task4_factor_overview(FactorSummary, tag, fp_fig_over);

        % Round split overview (optional)
        if ~isempty(FactorRoundSummary) && height(FactorRoundSummary)>0
            fp_csv_r = fullfile(fp_tbl_over, sprintf('task4_factor_round_summary_%s.csv', tag));
            writetable(FactorRoundSummary, fp_csv_r);
            % lightweight reuse: plot overviews separately per round
            try
                for rb = [1 2]
                    Trr = FactorRoundSummary(double(FactorRoundSummary.round)==rb,:);
                    if isempty(Trr) || height(Trr)==0, continue; end
                    fp_fig_over_r = fullfile(fp_fig_over, sprintf('round%d', rb));
                    mk(fp_fig_over_r);
                    plot_task4_factor_overview(Trr, sprintf('%s_round%d', tag, rb), fp_fig_over_r);
                end
            catch
            end
        end
    end
catch ME
    warning('analyze_core_lmm_suite: factor overview failed: %s', ME.message);
end

try
    if ~isempty(TrendSummary) && height(TrendSummary)>0
        fp_tbl_over_t = fullfile(fp_trend, 'tables', tag, 'overview');
        fp_fig_over_t = fullfile(fp_trend, 'figures', tag, 'overview');
        mk(fp_tbl_over_t); mk(fp_fig_over_t);
        fp_csv_t = fullfile(fp_tbl_over_t, sprintf('task4_trend_summary_%s.csv', tag));
        writetable(TrendSummary, fp_csv_t);
        plot_task4_trend_overview(TrendSummary, tag, fp_fig_over_t);
        if ~isempty(TrendRoundSummary) && height(TrendRoundSummary)>0
            fp_csv_tr = fullfile(fp_tbl_over_t, sprintf('task4_trend_round_summary_%s.csv', tag));
            writetable(TrendRoundSummary, fp_csv_tr);
            plot_task4_trend_round_overview(TrendRoundSummary, tag, fp_fig_over_t);
        end
    end
catch ME
    warning('analyze_core_lmm_suite: trend overview failed: %s', ME.message);
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
    C = to_table_compat(lme.Coefficients);
    % Holm adjustment within this fixed-effects table (exclude intercept if present)
    try
        C2 = C;
        if all(ismember({'Name','pValue'}, C2.Properties.VariableNames))
            nm = string(C2.Name);
            isInt = (nm=="(Intercept)") | (lower(nm)=="intercept");
            p = double(C2.pValue);
            padj = nan(size(p));
            padj(~isInt) = pipeline.holm_stepdown(p(~isInt));
            C2.p_holm = padj;
        else
            C2 = pipeline.add_holm_to_anova_table(C2);
        end
        C = C2;
    catch
    end
    writetable(C, fullfile(fp_tbl, sprintf('%s_fixed_effects.csv', stem)));
catch ME
    fprintf(2,'[WARN] task4 write_lme_tables fixed failed (%s): %s\n', stem, ME.message);
end
% anova table
try
    A = to_table_compat(anova(lme,'DFMethod','Satterthwaite'));
    try
        A = pipeline.add_holm_to_anova_table(A);
    catch
    end
    writetable(A, fullfile(fp_tbl, sprintf('%s_anova.csv', stem)));
catch ME
    fprintf(2,'[WARN] task4 write_lme_tables anova failed (%s): %s\n', stem, ME.message);
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
legend(ax,'Location','northeastoutside');
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

function [pWWR, pCx, pGrp, pWxC, pWxG, pCxG, p3] = extract_factor_pvals(lme1, lme2, lme3, did3)
pWWR = NaN; pCx = NaN; pGrp = NaN; pWxC = NaN; pWxG = NaN; pCxG = NaN; p3 = NaN;

% Main effects must come from Model1 (design contract)
try
    A1 = to_table_compat(anova(lme1,'DFMethod','Satterthwaite'));
    if all(ismember({'Term','pValue'}, A1.Properties.VariableNames))
        tt1 = string(A1.Term);
        pWWR = get_term_p(tt1, A1.pValue, "WWR");
        pCx  = get_term_p(tt1, A1.pValue, "Complexity");
        pGrp = get_term_p(tt1, A1.pValue, "Group");
    else
        fprintf(2,'[WARN] task4 extract_factor_pvals: Model1 ANOVA missing Term/pValue columns.\n');
    end
catch ME
    fprintf(2,'[WARN] task4 extract_factor_pvals A1 failed: %s\n', ME.message);
end

% Two-way interactions come from Model2
try
    A2 = to_table_compat(anova(lme2,'DFMethod','Satterthwaite'));
    if all(ismember({'Term','pValue'}, A2.Properties.VariableNames))
        tt2 = string(A2.Term);
        pWxC = get_interaction_p(tt2, A2.pValue, ["WWR","Complexity"]);
        pWxG = get_interaction_p(tt2, A2.pValue, ["WWR","Group"]);
        pCxG = get_interaction_p(tt2, A2.pValue, ["Complexity","Group"]);
    else
        fprintf(2,'[WARN] task4 extract_factor_pvals: Model2 ANOVA missing Term/pValue columns.\n');
    end
catch ME
    fprintf(2,'[WARN] task4 extract_factor_pvals A2 failed: %s\n', ME.message);
end

if did3 && ~isempty(lme3)
    try
        A3 = to_table_compat(anova(lme3,'DFMethod','Satterthwaite'));
        if all(ismember({'Term','pValue'}, A3.Properties.VariableNames))
            tt3 = string(A3.Term);
            i3 = find(contains(tt3,"WWR") & contains(tt3,"Complexity") & contains(tt3,"Group"),1,'first');
            if ~isempty(i3), p3 = double(A3.pValue(i3)); end
        end
    catch ME
        fprintf(2,'[WARN] task4 extract_factor_pvals A3 failed: %s\n', ME.message);
    end
end
end

function p = get_term_p(tt, pv, name)
p = NaN;
idx = find(tt==name,1,'first');
if ~isempty(idx), p = double(pv(idx)); end
end

function p = get_interaction_p(tt, pv, names)
p = NaN;
idx = find(contains(tt,names(1)) & contains(tt,names(2)) & contains(tt,":"),1,'first');
if ~isempty(idx), p = double(pv(idx)); end
end

function TR = collect_round_trend_stats(Tt, analysisName, dv)
TR = table();
if ~ismember('block_id', Tt.Properties.VariableNames)
    return;
end

for b = [1 2]
    idx = isfinite(Tt.block_id) & (double(Tt.block_id)==b) & ...
          isfinite(Tt.WWRc) & isfinite(Tt.WWRc2) & isfinite(Tt.EEG);
    Tb = Tt(idx,:);
    if height(Tb) < 18
        continue;
    end

    try
        lmeLb = fitlme(Tb, 'EEG ~ WWRc + Complexity + Group + (1|Subject)');
        lmeQb = fitlme(Tb, 'EEG ~ WWRc + WWRc2 + Complexity + Group + (1|Subject)');
        [bL, pL, bQ2, pQ2] = extract_trend_stats(lmeLb, lmeQb);
        isInvU = isfinite(bQ2) && isfinite(pQ2) && (bQ2 < 0) && (pQ2 < 0.05);
        rr = table(string(analysisName), string(dv), b, height(Tb), numel(unique(string(Tb.Subject))), ...
            bL, pL, bQ2, pQ2, isInvU, ...
            'VariableNames', {'analysis','metric','round','n_rows','n_subjects','beta_linear','p_linear','beta_quadratic','p_quadratic','is_inverted_u'});
        TR = [TR; rr]; %#ok<AGROW>
    catch
    end
end
end

function [bL, pL, bQ2, pQ2] = extract_trend_stats(lmeL, lmeQ)
bL = NaN; pL = NaN; bQ2 = NaN; pQ2 = NaN;
try
    CL = to_table_compat(lmeL.Coefficients); nL = string(CL.Name);
    iL = find(nL=="WWRc",1,'first');
    if ~isempty(iL)
        bL = double(CL.Estimate(iL));
        pL = double(CL.pValue(iL));
    end
    AL = to_table_compat(anova(lmeL,'DFMethod','Satterthwaite'));
    if all(ismember({'Term','pValue'}, AL.Properties.VariableNames))
        iLa = find(string(AL.Term)=="WWRc",1,'first');
        if ~isempty(iLa), pL = double(AL.pValue(iLa)); end
    end
catch ME
    fprintf(2,'[WARN] task4 extract_trend_stats linear failed: %s\n', ME.message);
end
try
    CQ = to_table_compat(lmeQ.Coefficients); nQ = string(CQ.Name);
    iQ = find(nQ=="WWRc2",1,'first');
    if ~isempty(iQ)
        bQ2 = double(CQ.Estimate(iQ));
        pQ2 = double(CQ.pValue(iQ));
    end
    AQ = to_table_compat(anova(lmeQ,'DFMethod','Satterthwaite'));
    if all(ismember({'Term','pValue'}, AQ.Properties.VariableNames))
        iQa = find(string(AQ.Term)=="WWRc2",1,'first');
        if ~isempty(iQa), pQ2 = double(AQ.pValue(iQa)); end
    end
catch ME
    fprintf(2,'[WARN] task4 extract_trend_stats quadratic failed: %s\n', ME.message);
end
end

function plot_task4_factor_overview(FS, tag, fp_fig_over)
analyses = unique(string(FS.analysis),'stable');
if isempty(analyses), return; end
for ai=1:numel(analyses)
    an = analyses(ai);
    T = FS(string(FS.analysis)==an,:);
    mets = unique(string(T.metric),'stable');
    cols = {'p_WWR','p_Complexity','p_Group','p_WWRxComplexity','p_WWRxGroup','p_ComplexityxGroup','p_threeway'};
    labels = {'WWR','Cx','Group','WWR×Cx','WWR×Group','Cx×Group','3-way'};
    Z = nan(numel(mets), numel(cols));
    for r=1:numel(mets)
        tr = T(string(T.metric)==mets(r),:);
        if isempty(tr), continue; end
        for c=1:numel(cols)
            if ismember(cols{c}, tr.Properties.VariableNames)
                Z(r,c) = double(tr.(cols{c})(1));
            end
        end
    end
    try
        nZero = sum(Z(:)==0, 'omitnan');
        if nZero > 0
            fprintf(2,'[WARN] task4 factor overview (%s/%s): detected %d zero p-values; display will clamp p at floor for color scaling.\n', tag, an, nZero);
        end
    catch
    end
    % Significance-tier heatmap for robust readability across MATLAB versions:
    % 0 = ns (p>=0.05), 1 = * (p<0.05), 2 = ** (p<0.01), 3 = *** (p<0.001)
    V = nan(size(Z));
    V(isfinite(Z)) = 0;
    V(isfinite(Z) & Z<0.05) = 1;
    V(isfinite(Z) & Z<0.01) = 2;
    V(isfinite(Z) & Z<0.001) = 3;

    set(0,'DefaultFigureVisible','off');
    fig = figure('Color','w','Position',[90 90 1180 420]);
    ax = axes(fig); hold(ax,'on');
    hImg = imagesc(ax, V);
    hImg.AlphaData = ~isnan(V);
    set(ax,'YDir','normal'); axis(ax,'tight');
    colormap(ax, interp1([0 0.5 1], [0.72 0.34 0.26; 0.96 0.96 0.96; 0.31 0.47 0.67], linspace(0,1,256)));
    cb = colorbar(ax); cb.Label.String = 'significance tier (0=ns,1=*,2=**,3=***)';
    caxis(ax, [0 3]);
    set(ax,'XTick',1:numel(cols),'XTickLabel',labels);
    set(ax,'YTick',1:numel(mets),'YTickLabel',cellstr(mets));
    xtickangle(ax,20);
    title(ax, sprintf('Task4 factor overview | %s [%s]', an, tag), 'Interpreter','none');

    for r=1:size(Z,1)
        for c=1:size(Z,2)
            if isnan(Z(r,c)), continue; end
            star=''; p=Z(r,c);
            if p<0.001, star='***'; elseif p<0.01, star='**'; elseif p<0.05, star='*'; end
            if p==0
                ptxt = sprintf('p<1e-300%s', star);
            else
                ptxt = sprintf('p=%.3g%s', p, star);
            end
            text(ax,c,r,ptxt,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'Color',[0.1 0.1 0.1]);
        end
    end

    fp = fullfile(fp_fig_over, pipeline.sanitize_filename(sprintf('task4_factor_overview_%s_%s.png', tag, an)));
    pipeline.export_figure_png(fig, fp, 300);
    try; close(fig); catch; end
end
end

function plot_task4_trend_round_overview(TRS, tag, fp_fig_over_t)
analyses = unique(string(TRS.analysis),'stable');
if isempty(analyses), return; end
for ai=1:numel(analyses)
    an = analyses(ai);
    T = TRS(string(TRS.analysis)==an,:);
    mets = unique(string(T.metric),'stable');
    Z = nan(2, numel(mets));
    P = nan(2, numel(mets));
    B = nan(2, numel(mets));
    for c=1:numel(mets)
        tr1 = T(string(T.metric)==mets(c) & double(T.round)==1,:);
        tr2 = T(string(T.metric)==mets(c) & double(T.round)==2,:);
        if ~isempty(tr1)
            B(1,c) = double(tr1.beta_quadratic(1));
            P(1,c) = double(tr1.p_quadratic(1));
            Z(1,c) = double(tr1.is_inverted_u(1));
        end
        if ~isempty(tr2)
            B(2,c) = double(tr2.beta_quadratic(1));
            P(2,c) = double(tr2.p_quadratic(1));
            Z(2,c) = double(tr2.is_inverted_u(1));
        end
    end

    set(0,'DefaultFigureVisible','off');
    fig = figure('Color','w','Position',[90 90 1180 400]);
    ax = axes(fig); hold(ax,'on');
    imagesc(ax, Z); set(ax,'YDir','normal'); axis(ax,'tight');
    colormap(ax, [0.94 0.94 0.94; 0.22 0.63 0.53]);
    cb = colorbar(ax); cb.Ticks = [0 1]; cb.TickLabels = {'not Inverted-U','Inverted-U'}; cb.Label.String = 'quadratic verdict';
    caxis(ax,[0 1]);
    set(ax,'XTick',1:numel(mets),'XTickLabel',cellstr(mets));
    set(ax,'YTick',[1 2],'YTickLabel',{'Round1 / Block1','Round2 / Block2'});
    xtickangle(ax,20);
    title(ax, sprintf('Task4 inverted-U by round | %s [%s]', an, tag), 'Interpreter','none');

    for r=1:2
        for c=1:numel(mets)
            if isnan(Z(r,c)), continue; end
            p=P(r,c); b=B(r,c); star='';
            if ~isnan(p)
                if p<0.001, star='***'; elseif p<0.01, star='**'; elseif p<0.05, star='*'; end
            end
            text(ax,c,r,sprintf('b2=%.3g\np=%.3g%s',b,p,star), ...
                'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'Color',[0.1 0.1 0.1]);
        end
    end

    fp = fullfile(fp_fig_over_t, pipeline.sanitize_filename(sprintf('task4_invertedu_round_overview_%s_%s.png', tag, an)));
    pipeline.export_figure_png(fig, fp, 300);
    try; close(fig); catch; end
end
end

function plot_task4_trend_overview(TS, tag, fp_fig_over_t)
analyses = unique(string(TS.analysis),'stable');
if isempty(analyses), return; end
for ai=1:numel(analyses)
    an = analyses(ai);
    T = TS(string(TS.analysis)==an,:);
    mets = unique(string(T.metric),'stable');
    Z = nan(2, numel(mets));
    P = nan(2, numel(mets));
    for c=1:numel(mets)
        tr = T(string(T.metric)==mets(c),:);
        if isempty(tr), continue; end
        % With imagesc + YDir='normal': row1 is bottom, row2 is top.
        % Force display order: top=linear, bottom=quadratic.
        Z(1,c) = double(tr.beta_quadratic(1));
        Z(2,c) = double(tr.beta_linear(1));
        P(1,c) = double(tr.p_quadratic(1));
        P(2,c) = double(tr.p_linear(1));
    end

    set(0,'DefaultFigureVisible','off');
    fig = figure('Color','w','Position',[90 90 1180 400]);
    ax = axes(fig); hold(ax,'on');
    imagesc(ax, Z); set(ax,'YDir','normal'); axis(ax,'tight');
    colormap(ax, interp1([0 0.5 1], [0.72 0.34 0.26; 0.96 0.96 0.96; 0.31 0.47 0.67], linspace(0,1,256)));
    cb = colorbar(ax); cb.Label.String = 'beta';
    mx = max(abs(Z(:)),[],'omitnan'); if isempty(mx)||~isfinite(mx)||mx==0, mx=0.01; end
    caxis(ax,[-mx mx]);
    set(ax,'XTick',1:numel(mets),'XTickLabel',cellstr(mets));
    set(ax,'YTick',[1 2],'YTickLabel',{'quadratic (WWRc2) [bottom]','linear (WWRc) [top]'});
    xtickangle(ax,20);
    title(ax, sprintf('Task4 trend overview | %s [%s]', an, tag), 'Interpreter','none');

    for r=1:2
        for c=1:numel(mets)
            if isnan(Z(r,c)), continue; end
            p=P(r,c); star='';
            if ~isnan(p)
                if p<0.001, star='***'; elseif p<0.01, star='**'; elseif p<0.05, star='*'; end
            end
            text(ax,c,r,sprintf('β=%.3g\np=%.3g%s',Z(r,c),p,star),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
        end
    end

    fp = fullfile(fp_fig_over_t, pipeline.sanitize_filename(sprintf('task4_trend_overview_%s_%s.png', tag, an)));
    pipeline.export_figure_png(fig, fp, 300);
    try; close(fig); catch; end
end
end

function Ttbl = to_table_compat(X)
if istable(X)
    Ttbl = X;
    return;
end
try
    if isa(X,'dataset')
        Ttbl = dataset2table(X);
        return;
    end
catch
end
try
    Ttbl = struct2table(X);
    return;
catch
end
error('Unsupported output type for table export: %s', class(X));
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
