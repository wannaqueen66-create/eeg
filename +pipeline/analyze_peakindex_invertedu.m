function out = analyze_peakindex_invertedu(AllScene, fp_sum, cfg, tag)
%ANALYZE_PEAKINDEX_INVERTEDU Task5: Inverted-U PeakIndex test.
%
% For each metric (default: O_theta, F_theta), and each grouping analysis
% (ExperienceGroup / SportFreqGroup), compute subject-level:
%   PeakIndex = EEG(45) - mean(EEG(15), EEG(75))
% within each Complexity level.
%
% Then fit:
%   PeakIndex ~ Complexity*Group + (1|Subject)
%
% Outputs under:
%   <summary>/analysis-2/task5_peakindex_invertedU/
%     tables/<tag>/{experience|sportfreq}/...
%     reports/<tag>/{experience|sportfreq}/...
%     figures/<tag>/{experience|sportfreq}/...

if nargin < 4 || isempty(tag)
    tag = 'raw';
end
if nargin < 3
    cfg = struct();
end

out = struct();

if exist('fitlme','file') ~= 2
    warning('analyze_peakindex_invertedu: fitlme not found (Stats toolbox missing). Skipping.');
    return;
end

req = {'subject_id','WWR','Complexity'};
for i=1:numel(req)
    if ~ismember(req{i}, AllScene.Properties.VariableNames)
        warning('analyze_peakindex_invertedu: missing required column %s. Skipping.', req{i});
        return;
    end
end

metrics = ["O_theta","F_theta"];
try
    if isfield(cfg,'task5_metrics') && ~isempty(cfg.task5_metrics)
        metrics = string(cfg.task5_metrics);
    end
catch
end

analyses = { ...
    struct('name','experience','gcol','ExperienceGroup'), ...
    struct('name','sportfreq','gcol','SportFreqGroup') ...
};

fp_root = fullfile(fp_sum, 'analysis-2', 'task5_peakindex_invertedU');

T0 = AllScene;
T0.subject_id = string(T0.subject_id);
T0.WWR = normalize_wwr(T0.WWR);
T0.Complexity = normalize_complexity(T0.Complexity);

for ai = 1:numel(analyses)
    A = analyses{ai};
    gcol = string(A.gcol);
    if ~ismember(gcol, T0.Properties.VariableNames)
        continue;
    end

    for mi = 1:numel(metrics)
        dv = metrics(mi);
        if ~ismember(dv, T0.Properties.VariableNames)
            continue;
        end

        T = T0(:, {'subject_id','WWR','Complexity',gcol,dv});
        T.Properties.VariableNames = {'subject_id','WWR','Complexity','GroupRaw','EEGraw'};
        T.EEG = double(T.EEGraw);
        T.Group = normalize_high_low(T.GroupRaw);
        T = T(~isnan(T.EEG), :);
        T = T(strlength(string(T.WWR))>0 & strlength(string(T.Complexity))>0 & strlength(string(T.Group))>0, :);
        T = T(ismember(string(T.WWR), ["15","45","75"]), :);
        if height(T) < 30
            continue;
        end

        T.Subject = categorical(string(T.subject_id));
        T.WWR = categorical(string(T.WWR), {'15','45','75'});
        T.Complexity = categorical(string(T.Complexity), {'ComplexityLow','ComplexityHigh'});
        T.Group = categorical(string(T.Group), {'Low','High'});

        % subject-level mean EEG per WWR within (Subject, Complexity, Group)
        [G, sid, grp, cx, wwr] = findgroups(T.Subject, T.Group, T.Complexity, T.WWR);
        mu = splitapply(@(x) mean(x,'omitnan'), T.EEG, G);
        Awwr = table(sid, grp, cx, string(wwr), mu, 'VariableNames', {'Subject','Group','Complexity','WWR','EEG'});

        % Wide by WWR level
        K = unique(Awwr(:, {'Subject','Group','Complexity'}), 'rows');
        for lv = ["15","45","75"]
            Z = Awwr(Awwr.WWR==lv, {'Subject','Group','Complexity','EEG'});
            Z.Properties.VariableNames{'EEG'} = ['EEG_' char(lv)];
            K = outerjoin(K, Z, 'Keys', {'Subject','Group','Complexity'}, 'MergeKeys', true, 'Type','left');
        end

        if ~ismember('EEG_45', K.Properties.VariableNames)
            continue;
        end
        if ~ismember('EEG_15', K.Properties.VariableNames)
            K.EEG_15 = nan(height(K),1);
        end
        if ~ismember('EEG_75', K.Properties.VariableNames)
            K.EEG_75 = nan(height(K),1);
        end

        K.PeakIndex = K.EEG_45 - mean([K.EEG_15 K.EEG_75], 2, 'omitnan');
        K = K(~isnan(K.PeakIndex), :);
        if height(K) < 12
            continue;
        end

        % Output dirs
        fp_tbl = fullfile(fp_root, 'tables', tag, A.name);
        fp_rep = fullfile(fp_root, 'reports', tag, A.name);
        fp_fig = fullfile(fp_root, 'figures', tag, A.name);
        if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
        if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
        if ~exist(fp_fig,'dir'); mkdir(fp_fig); end

        % Save subject-level PeakIndex
        fp_subject = fullfile(fp_tbl, sprintf('peakindex_subjectlevel_%s_%s.csv', dv, tag));
        writetable(K, fp_subject);

        % LMM
        try
            lme = fitlme(K, 'PeakIndex ~ Complexity*Group + (1|Subject)');
        catch ME
            warning('analyze_peakindex_invertedu: fitlme failed (%s/%s/%s): %s', tag, A.name, dv, ME.message);
            continue;
        end

        % export tables
        try
            writetable(lme.Coefficients, fullfile(fp_tbl, sprintf('peakindex_lmm_fixed_effects_%s_%s.csv', dv, tag)));
        catch
        end
        try
            Aov = anova(lme,'DFMethod','Satterthwaite');
            writetable(Aov, fullfile(fp_tbl, sprintf('peakindex_lmm_anova_%s_%s.csv', dv, tag)));
        catch
        end

        % report
        try
            fp_md = fullfile(fp_rep, sprintf('peakindex_report_%s_%s.md', dv, tag));
            write_report(fp_md, dv, tag, A.name, K, lme);
        catch
        end

        % figure
        try
            fp_png = fullfile(fp_fig, pipeline.sanitize_filename(sprintf('peakindex_%s_%s_%s.png', tag, A.name, dv)));
            plot_peakindex_figure(K, dv, tag, A.name, fp_png);
        catch
        end
    end
end

end

function plot_peakindex_figure(K, dv, tag, analysisName, fp_png)
set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[90 90 900 380]);

levelsCx = categories(K.Complexity);
if isempty(levelsCx)
    levelsCx = {'ComplexityLow','ComplexityHigh'};
end

for i=1:min(2,numel(levelsCx))
    cx = string(levelsCx{i});
    ax = subplot(1,2,i); hold(ax,'on');

    D = K(string(K.Complexity)==cx, :);
    groups = ["Low","High"];
    cols = [0.2 0.5 0.9; 0.9 0.4 0.2];

    for gi=1:2
        g = groups(gi);
        y = D.PeakIndex(string(D.Group)==g);
        if isempty(y); continue; end
        x = gi + (rand(size(y))-0.5)*0.14;
        scatter(ax, x, y, 24, cols(gi,:), 'filled', 'MarkerFaceAlpha',0.6, 'HandleVisibility','off');
        mu = mean(y,'omitnan');
        se = std(y,'omitnan')/sqrt(numel(y));
        errorbar(ax, gi, mu, se, 'o', 'Color', cols(gi,:), 'MarkerFaceColor', cols(gi,:), 'LineWidth',1.4, 'DisplayName', char(g));
    end

    xlim(ax,[0.5 2.5]);
    set(ax,'XTick',1:2,'XTickLabel',{'Low','High'});
    ylabel(ax,'PeakIndex = EEG45 - mean(EEG15, EEG75)');
    title(ax, char(cx), 'Interpreter','none');
    grid(ax,'on');
end

sgtitle(sprintf('Task5 PeakIndex | %s | %s [%s]', analysisName, dv, tag), 'Interpreter','none');
pipeline.export_figure_png(fig, fp_png, 300);
try; close(fig); catch; end
end

function write_report(fp_md, dv, tag, analysisName, K, lme)
lines = strings(0,1);
lines(end+1) = sprintf('# Task5 PeakIndex (Inverted-U) | %s | %s | %s', analysisName, dv, tag);
lines(end+1) = '';
lines(end+1) = 'Model: `PeakIndex ~ Complexity*Group + (1|Subject)`';
lines(end+1) = '';

% quick means
try
    lines(end+1) = '## Raw means by Complexity × Group';
    lines(end+1) = 'Complexity | Group | n | mean | sem';
    lines(end+1) = '---|---:|---:|---:|---:';
    [G, cx, gr] = findgroups(string(K.Complexity), string(K.Group));
    n = splitapply(@numel, K.PeakIndex, G);
    mu = splitapply(@(x) mean(x,'omitnan'), K.PeakIndex, G);
    se = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), K.PeakIndex, G);
    for i=1:numel(n)
        lines(end+1) = sprintf('%s|%s|%d|%.6g|%.6g', cx(i), gr(i), n(i), mu(i), se(i));
    end
    lines(end+1) = '';
catch
end

lines(end+1) = '## LMM outputs';
lines(end+1) = '- See tables: `peakindex_lmm_fixed_effects_*.csv` and `peakindex_lmm_anova_*.csv`.';
lines(end+1) = '- Interpret Complexity and Group main effects first, then Complexity×Group.';

try
    C = lme.Coefficients;
    rn = string(C.Name);
    ix = find(contains(rn,'Complexity') & contains(rn,'Group'),1);
    if ~isempty(ix)
        lines(end+1) = sprintf('- Complexity×Group p = %.4g', C.pValue(ix));
    end
catch
end

fid = fopen(fp_md,'w');
for i=1:numel(lines)
    fprintf(fid, '%s\n', lines(i));
end
fclose(fid);
end

function w = normalize_wwr(x)
w = string(x);
w = strtrim(w);
for i=1:numel(w)
    tok = regexp(char(w(i)), '(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        w(i) = string(str2double(tok{1}));
    end
end
ok = ismember(w,["15","45","75"]);
w(~ok) = "";
end

function c = normalize_complexity(x)
c = string(x);
c = strtrim(c);
cl = lower(c);
out = repmat("", numel(c), 1);
out(ismember(cl,["low","0","c0","complexitylow"])) = "ComplexityLow";
out(ismember(cl,["high","1","c1","complexityhigh"])) = "ComplexityHigh";
out(c=="ComplexityLow") = "ComplexityLow";
out(c=="ComplexityHigh") = "ComplexityHigh";
isNum = ~isnan(str2double(c));
out(isNum & str2double(c)==0) = "ComplexityLow";
out(isNum & str2double(c)==1) = "ComplexityHigh";
c = out;
end

function g = normalize_high_low(x)
s = string(x);
s = strtrim(s);
sl = lower(s);
g = repmat("", numel(s), 1);
g(ismember(sl,["high","1","高","h"])) = "High";
g(ismember(sl,["low","0","低","l"])) = "Low";
mask = (s=="High" | s=="Low");
g(mask) = s(mask);
end
