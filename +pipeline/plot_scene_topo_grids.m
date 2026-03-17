function plot_scene_topo_grids(mode, EEG, T, bands, totalBand, fs, wlen, nover, nfft, fp_fig, fp_csv, fp_qc, base, designMap, fp_batch, cfg)
%PLOT_SCENE_TOPO_GRIDS Render scene-level topo grids arranged by WWR x Complexity.
%
% mode="subject": build per-subject scene topo grids (2x3 per block, 2 blocks total)
%   and export long-form channel values for later group aggregation.
%
% mode="group": read aggregated per-subject scene topo exports and render
%   group-mean scene topo grids under the batch figure tree.
%
% Expected scene ordering per block:
%   row1 = ComplexityLow  + WWR 15 / 45 / 75
%   row2 = ComplexityHigh + WWR 15 / 45 / 75
%
% For the 12 scenes, figures are produced when data exist:
%   block1    -> scenes 1..6
%   block2    -> scenes 7..12
%   allblocks -> collapse block1+block2 into one 2x3 grid by WWR x Complexity

if nargin < 1 || isempty(mode)
    mode = 'subject';
end
mode = lower(strtrim(string(mode)));

switch mode
    case "subject"
        plot_subject_scene_topos(EEG, T, bands, totalBand, fs, wlen, nover, nfft, fp_fig, fp_csv, fp_qc, base, designMap);
    case "group"
        plot_group_scene_topos(EEG, bands, fp_fig, fp_batch, cfg);
    otherwise
        error('plot_scene_topo_grids: unsupported mode %s', mode);
end
end

function plot_subject_scene_topos(EEG, T, bands, totalBand, fs, wlen, nover, nfft, fp_fig, fp_csv, fp_qc, base, designMap)
if ~isfield(EEG,'chanlocs') || numel(EEG.chanlocs)~=EEG.nbchan
    warning('plot_scene_topo_grids(subject): no valid chanlocs found.');
    return;
end
if isempty(T) || ~ismember('scene_id', T.Properties.VariableNames)
    warning('plot_scene_topo_grids(subject): missing scene_id.');
    return;
end
if nargin < 13 || isempty(designMap)
    designMap = table();
end

try
    pipeline.write_chanlocs_snapshot(fp_qc, base, EEG.chanlocs);
catch
end

S = table((1:height(T))', double(T.scene_id), string(T.cond), 'VariableNames', {'seg_idx','scene_id','cond'});
S = S(lower(strtrim(S.cond))=="view", :);
if isempty(S)
    warning('plot_scene_topo_grids(subject): no view segments.');
    return;
end
try
    S = pipeline.attach_design(S, base, designMap);
catch
end
if ~ismember('WWR', S.Properties.VariableNames) || ~ismember('Complexity', S.Properties.VariableNames)
    warning('plot_scene_topo_grids(subject): missing WWR/Complexity after design attachment.');
    return;
end

S.WWRn = normalize_wwr_local(S.WWR);
S.CX = normalize_complexity_local(S.Complexity);
S.Block = nan(height(S),1);
if ismember('block_id', T.Properties.VariableNames)
    try
        S.Block = double(T.block_id(S.seg_idx));
    catch
    end
end
if all(~isfinite(S.Block))
    S.Block = infer_block_from_scene_id(S.scene_id);
end
S = S(isfinite(S.scene_id) & isfinite(S.Block) & strlength(S.WWRn)>0 & strlength(S.CX)>0, :);
if isempty(S)
    warning('plot_scene_topo_grids(subject): no usable scene/design rows after filtering.');
    return;
end

bandNames = {'theta','alpha','beta'};
bandRanges = {bands.theta, bands.alpha, bands.beta};
allRows = {};

for bi = 1:3
    bandName = string(bandNames{bi});
    bandRange = bandRanges{bi};

    [sceneTopo, meta] = compute_scene_topos(EEG, T, S, bandRange, totalBand, fs, wlen, nover, nfft);
    if isempty(sceneTopo)
        continue;
    end

    for ii = 1:numel(sceneTopo)
        allRows(end+1,:) = {string(base), double(meta.scene_id(ii)), double(meta.Block(ii)), string(meta.WWRn(ii)), string(meta.CX(ii)), bandName, sceneTopo{ii}}; %#ok<AGROW>
    end

    for blk = 1:2
        use = find(meta.Block==blk);
        if isempty(use), continue; end
        fig = figure('Color','w','Position',[100 100 1350 760]);
        sgtitle(sprintf('%s scene topoplots | %s | block%d', upper(char(bandName)), base, blk), 'Interpreter','none', 'FontWeight','normal');
        orderTbl = make_scene_grid_order(meta(use,:));
        try
            writetable(add_layout_indices(orderTbl, blk), fullfile(fp_csv, sprintf('%s_scene_topogrid_layout_block%d.csv', base, blk)));
        catch
        end
        for k = 1:height(orderTbl)
            ax = subplot(2, 3, k); %#ok<LAXES>
            idx = find(meta.scene_id==orderTbl.scene_id(k) & meta.Block==blk & meta.WWRn==orderTbl.WWRn(k) & meta.CX==orderTbl.CX(k), 1, 'first');
            if isempty(idx)
                axis off;
                title(sprintf('missing | WWR%s %s', char(orderTbl.WWRn(k)), char(orderTbl.CX(k))), 'Interpreter','none');
                continue;
            end
            topoplot(sceneTopo{idx}, EEG.chanlocs, 'electrodes','labels');
            title(sprintf('scene%02d | WWR%s | %s', orderTbl.scene_id(k), char(orderTbl.WWRn(k)), short_cx(orderTbl.CX(k))), 'Interpreter','none', 'FontSize',10);
            colorbar;
        end
        fn = pipeline.sanitize_filename(sprintf('%s_scene_topogrid_%s_block%d.png', base, bandName, blk));
        saveas(fig, fullfile(fp_fig, fn));
        try; close(fig); catch; end
    end

    [comboTopo, comboMeta] = collapse_scene_topos_across_blocks(sceneTopo, meta);
    if ~isempty(comboTopo)
        fig = figure('Color','w','Position',[100 100 1350 760]);
        sgtitle(sprintf('%s scene topoplots | %s | allblocks', upper(char(bandName)), base), 'Interpreter','none', 'FontWeight','normal');
        orderTbl = make_scene_grid_order(comboMeta);
        try
            writetable(add_layout_indices(orderTbl, NaN), fullfile(fp_csv, sprintf('%s_scene_topogrid_layout_allblocks.csv', base)));
        catch
        end
        for k = 1:height(orderTbl)
            ax = subplot(2, 3, k); %#ok<LAXES>
            idx = find(comboMeta.WWRn==orderTbl.WWRn(k) & comboMeta.CX==orderTbl.CX(k), 1, 'first');
            if isempty(idx)
                axis off;
                title(sprintf('missing | WWR%s %s', char(orderTbl.WWRn(k)), char(orderTbl.CX(k))), 'Interpreter','none');
                continue;
            end
            topoplot(comboTopo{idx}, EEG.chanlocs, 'electrodes','labels');
            title(sprintf('combined | WWR%s | %s', char(orderTbl.WWRn(k)), short_cx(orderTbl.CX(k))), 'Interpreter','none', 'FontSize',10);
            colorbar;
        end
        fn = pipeline.sanitize_filename(sprintf('%s_scene_topogrid_%s_allblocks.png', base, bandName));
        saveas(fig, fullfile(fp_fig, fn));
        try; close(fig); catch; end
    end
end

if ~isempty(allRows)
    write_scene_topo_long(fullfile(fp_csv, sprintf('%s_scene_topo_long.csv', base)), allRows, EEG.nbchan);
end
end

function plot_group_scene_topos(EEG, bands, fp_fig, fp_batch, cfg)
if nargin < 6
    cfg = struct();
end
if exist('topoplot','file') ~= 2
    warning('plot_scene_topo_grids(group): topoplot not found on path.');
    return;
end
if nargin < 4 || isempty(fp_batch) || ~exist(fp_batch,'dir')
    warning('plot_scene_topo_grids(group): missing batch dir.');
    return;
end

chanlocs = [];
subjRoot = '';
try
    subjRoot = resolve_subjects_root_from_batch(fp_batch);
    d = dir(fullfile(subjRoot, '*', 'qc', '*_chanlocs.mat'));
    if ~isempty(d)
        S = load(fullfile(d(1).folder, d(1).name));
        if isfield(S,'chanlocs')
            chanlocs = S.chanlocs;
        end
    end
catch
end
if isempty(chanlocs) && isstruct(EEG) && isfield(EEG,'chanlocs')
    chanlocs = EEG.chanlocs;
end
if isempty(chanlocs)
    warning('plot_scene_topo_grids(group): no chanlocs available. resolved_subjects_root=%s', string(subjRoot));
    return;
end

if strlength(string(subjRoot))==0 || ~exist(subjRoot,'dir')
    subjRoot = resolve_subjects_root_from_batch(fp_batch);
end
D = dir(fullfile(subjRoot, '*', 'tables', '*_scene_topo_long.csv'));
if isempty(D)
    D = dir(fullfile(subjRoot, '*', 'csv', '*_scene_topo_long.csv'));
end
if isempty(D)
    warning('plot_scene_topo_grids(group): no per-subject scene_topo_long files found.');
    return;
end

Tall = table();
for i = 1:numel(D)
    try
        T = readtable(fullfile(D(i).folder, D(i).name), 'TextType','string');
        Tall = [Tall; T]; %#ok<AGROW>
    catch
    end
end
if isempty(Tall)
    return;
end

if ismember('subject_id', Tall.Properties.VariableNames)
    Tall.sid_key = canonical_subject_id_local(Tall.subject_id);
else
    warning('plot_scene_topo_grids(group): scene topo table missing subject_id.');
    return;
end

if isfield(cfg,'qc_include_subjects') && ~isempty(cfg.qc_include_subjects)
    try
        inc = canonical_subject_id_local(cfg.qc_include_subjects);
        Tall = Tall(ismember(Tall.sid_key, inc), :);
    catch
    end
end
if isempty(Tall)
    return;
end

if ~ismember('ExperienceGroup', Tall.Properties.VariableNames)
    Tall.ExperienceGroup = repmat("", height(Tall), 1);
else
    Tall.ExperienceGroup = normalize_high_low_local(Tall.ExperienceGroup);
end
try
    needAttach = all(strlength(strtrim(string(Tall.ExperienceGroup)))==0);
    if needAttach
        fp_tbl_raw = fp_batch;
        if exist('pipeline.get_table_dir','file')==2
            fp_tbl_raw = pipeline.get_table_dir(fp_batch, cfg, 'merged_raw');
        end
        f_scene = fullfile(fp_tbl_raw, 'all_subjects_scene_level.csv');
        if ~exist(f_scene,'file')
            f_scene = fullfile(fp_batch, 'all_subjects_scene_level.csv');
        end
        if exist(f_scene,'file')
            M = readtable(f_scene, 'TextType','string');
            if ismember('subject_id', M.Properties.VariableNames)
                sidm = canonical_subject_id_local(M.subject_id);
                [~, ia] = unique(sidm, 'stable');
                M = M(ia,:);
                sidm = sidm(ia);
                [tf, loc] = ismember(Tall.sid_key, sidm);
                src = repmat("", height(Tall), 1);
                if ismember('ExperienceGroup', M.Properties.VariableNames)
                    src(tf) = string(M.ExperienceGroup(loc(tf)));
                elseif ismember('Experience', M.Properties.VariableNames)
                    src(tf) = string(M.Experience(loc(tf)));
                end
                src = normalize_high_low_local(src);
                ex = Tall.ExperienceGroup;
                fillMask = tf & strlength(strtrim(src))>0;
                ex(fillMask) = src(fillMask);
                Tall.ExperienceGroup = ex;
            end
        end
    end
catch
end

bandNames = {'theta','alpha','beta'};
for bi = 1:3
    bandName = string(bandNames{bi});
    Tb = Tall(lower(strtrim(string(Tall.band)))==bandName, :);
    if isempty(Tb), continue; end
    meta = unique(Tb(:, {'scene_id','block_id','WWR','Complexity'}), 'rows');
    meta.WWRn = normalize_wwr_local(meta.WWR);
    meta.CX = normalize_complexity_local(meta.Complexity);
    meta = meta(strlength(meta.WWRn)>0 & strlength(meta.CX)>0, :);
    if isempty(meta), continue; end

    Tb.ExperienceGroup = normalize_high_low_local(Tb.ExperienceGroup);
    try
        fprintf('[scene topo group] band=%s | overall_rows=%d | exp_low_rows=%d | exp_high_rows=%d\n', ...
            char(bandName), height(Tb), sum(string(Tb.ExperienceGroup)=="Low"), sum(string(Tb.ExperienceGroup)=="High"));
    catch
    end

    groupDefs = {
        struct('name','overall','mode','mean','mask',true(height(Tb),1)), ...
        struct('name','experience_low','mode','mean','mask',string(Tb.ExperienceGroup)=="Low"), ...
        struct('name','experience_high','mode','mean','mask',string(Tb.ExperienceGroup)=="High"), ...
        struct('name','experience_highminuslow','mode','diff','maskH',string(Tb.ExperienceGroup)=="High",'maskL',string(Tb.ExperienceGroup)=="Low") ...
    };

    for gd = 1:numel(groupDefs)
        Gd = groupDefs{gd};
        for blk = 1:2
            useMeta = meta(double(meta.block_id)==blk,:);
            if isempty(useMeta), continue; end
            fig = figure('Color','w','Position',[100 100 1350 760]);
            sgtitle(sprintf('%s scene topoplots | %s | block%d', strrep(Gd.name,'_',' '), upper(char(bandName)), blk), 'Interpreter','none', 'FontWeight','normal');
            orderTbl = make_scene_grid_order(useMeta);
            try
                writetable(add_layout_indices(orderTbl, blk), fullfile(fp_fig, sprintf('%s_scene_topogrid_layout_block%d.csv', Gd.name, blk)));
            catch
            end
            for k = 1:height(orderTbl)
                ax = subplot(2, 3, k); %#ok<LAXES>
                baseMask = double(Tb.scene_id)==double(orderTbl.scene_id(k)) & double(Tb.block_id)==blk & normalize_wwr_local(Tb.WWR)==orderTbl.WWRn(k) & normalize_complexity_local(Tb.Complexity)==orderTbl.CX(k);
                vec = nan(numel(chanlocs),1);
                if strcmp(Gd.mode,'mean')
                    idxRows = baseMask & Gd.mask;
                    if any(idxRows)
                        G = groupsummary(Tb(idxRows,:), 'chan_idx', 'mean', 'value');
                        idxChan = double(G.chan_idx);
                        good = idxChan>=1 & idxChan<=numel(chanlocs);
                        vec(idxChan(good)) = double(G.mean_value(good));
                    end
                else
                    idxH = baseMask & Gd.maskH;
                    idxL = baseMask & Gd.maskL;
                    if any(idxH) && any(idxL)
                        GH = groupsummary(Tb(idxH,:), 'chan_idx', 'mean', 'value');
                        GL = groupsummary(Tb(idxL,:), 'chan_idx', 'mean', 'value');
                        vh = nan(numel(chanlocs),1); vl = nan(numel(chanlocs),1);
                        idxChan = double(GH.chan_idx); good = idxChan>=1 & idxChan<=numel(chanlocs); vh(idxChan(good)) = double(GH.mean_value(good));
                        idxChan = double(GL.chan_idx); good = idxChan>=1 & idxChan<=numel(chanlocs); vl(idxChan(good)) = double(GL.mean_value(good));
                        vec = vh - vl;
                    end
                end
                if all(~isfinite(vec))
                    try
                        fprintf(2, '[scene topo group] missing panel: group=%s band=%s block=%d WWR=%s CX=%s\n', ...
                            Gd.name, char(bandName), blk, char(orderTbl.WWRn(k)), char(orderTbl.CX(k)));
                    catch
                    end
                    axis off;
                    title(sprintf('missing | WWR%s %s', char(orderTbl.WWRn(k)), char(orderTbl.CX(k))), 'Interpreter','none');
                    continue;
                end
                topoplot(vec, chanlocs, 'electrodes','labels');
                title(sprintf('scene%02d | WWR%s | %s', orderTbl.scene_id(k), char(orderTbl.WWRn(k)), short_cx(orderTbl.CX(k))), 'Interpreter','none', 'FontSize',10);
                colorbar;
            end
            fn = pipeline.sanitize_filename(sprintf('%s_scene_topogrid_%s_block%d.png', Gd.name, bandName, blk));
            saveas(fig, fullfile(fp_fig, fn));
            try; close(fig); catch; end
        end

        comboMeta = unique(Tb(:, {'WWR','Complexity'}), 'rows');
        comboMeta.WWRn = normalize_wwr_local(comboMeta.WWR);
        comboMeta.CX = normalize_complexity_local(comboMeta.Complexity);
        comboMeta = comboMeta(strlength(comboMeta.WWRn)>0 & strlength(comboMeta.CX)>0, :);
        if ~isempty(comboMeta)
            fig = figure('Color','w','Position',[100 100 1350 760]);
            sgtitle(sprintf('%s scene topoplots | %s | allblocks', strrep(Gd.name,'_',' '), upper(char(bandName))), 'Interpreter','none', 'FontWeight','normal');
            orderTbl = make_scene_grid_order(comboMeta);
            try
                writetable(add_layout_indices(orderTbl, NaN), fullfile(fp_fig, sprintf('%s_scene_topogrid_layout_allblocks.csv', Gd.name)));
            catch
            end
            for k = 1:height(orderTbl)
                ax = subplot(2, 3, k); %#ok<LAXES>
                baseMask = normalize_wwr_local(Tb.WWR)==orderTbl.WWRn(k) & normalize_complexity_local(Tb.Complexity)==orderTbl.CX(k);
                vec = nan(numel(chanlocs),1);
                if strcmp(Gd.mode,'mean')
                    idxRows = baseMask & Gd.mask;
                    if any(idxRows)
                        G = groupsummary(Tb(idxRows,:), 'chan_idx', 'mean', 'value');
                        idxChan = double(G.chan_idx);
                        good = idxChan>=1 & idxChan<=numel(chanlocs);
                        vec(idxChan(good)) = double(G.mean_value(good));
                    end
                else
                    idxH = baseMask & Gd.maskH;
                    idxL = baseMask & Gd.maskL;
                    if any(idxH) && any(idxL)
                        GH = groupsummary(Tb(idxH,:), 'chan_idx', 'mean', 'value');
                        GL = groupsummary(Tb(idxL,:), 'chan_idx', 'mean', 'value');
                        vh = nan(numel(chanlocs),1); vl = nan(numel(chanlocs),1);
                        idxChan = double(GH.chan_idx); good = idxChan>=1 & idxChan<=numel(chanlocs); vh(idxChan(good)) = double(GH.mean_value(good));
                        idxChan = double(GL.chan_idx); good = idxChan>=1 & idxChan<=numel(chanlocs); vl(idxChan(good)) = double(GL.mean_value(good));
                        vec = vh - vl;
                    end
                end
                if all(~isfinite(vec))
                    try
                        fprintf(2, '[scene topo group] missing panel: group=%s band=%s allblocks WWR=%s CX=%s\n', ...
                            Gd.name, char(bandName), char(orderTbl.WWRn(k)), char(orderTbl.CX(k)));
                    catch
                    end
                    axis off;
                    title(sprintf('missing | WWR%s %s', char(orderTbl.WWRn(k)), char(orderTbl.CX(k))), 'Interpreter','none');
                    continue;
                end
                topoplot(vec, chanlocs, 'electrodes','labels');
                title(sprintf('combined | WWR%s | %s', char(orderTbl.WWRn(k)), short_cx(orderTbl.CX(k))), 'Interpreter','none', 'FontSize',10);
                colorbar;
            end
            fn = pipeline.sanitize_filename(sprintf('%s_scene_topogrid_%s_allblocks.png', Gd.name, bandName));
            saveas(fig, fullfile(fp_fig, fn));
            try; close(fig); catch; end
        end
    end
end
end

function [sceneTopo, meta] = compute_scene_topos(EEG, T, S, bandRange, totalBand, fs, wlen, nover, nfft)
meta = unique(S(:, {'scene_id','Block','WWRn','CX'}), 'rows');
meta = sortrows(meta, {'Block','scene_id'});
sceneTopo = cell(height(meta),1);
for i = 1:height(meta)
    use = S.scene_id==meta.scene_id(i) & S.Block==meta.Block(i) & S.WWRn==meta.WWRn(i) & S.CX==meta.CX(i);
    idxSeg = S.seg_idx(use);
    if isempty(idxSeg), continue; end
    sceneTopo{i} = mean_rel_topo_local(EEG, T, idxSeg, bandRange, totalBand, fs, wlen, nover, nfft);
end
end

function rB = mean_rel_topo_local(EEG, T, idxSeg, bandRange, totalBand, fs, wlen, nover, nfft)
rB = zeros(EEG.nbchan,1);
cnt = 0;
for kk = idxSeg(:)'
    seg = double(EEG.data(:, T.s0(kk):T.s1(kk)));
    [P,F] = pwelch(seg', wlen, nover, nfft, fs);
    for ci = 1:EEG.nbchan
        [rel, ~] = band_power_local(P(:,ci), F, bandRange, totalBand);
        rB(ci) = rB(ci) + rel;
    end
    cnt = cnt + 1;
end
rB = rB / max(cnt,1);
end

function [bp, ib] = band_power_local(Pxx, F, bandRange, totalBand)
ib = F >= bandRange(1) & F <= bandRange(2);
it = F >= totalBand(1) & F <= totalBand(2);
bpBand = trapz(F(ib), Pxx(ib));
bpTot  = trapz(F(it), Pxx(it));
bp = bpBand / max(bpTot, eps);
end

function write_scene_topo_long(fp, rows, nbchan)
out = {};
for i = 1:size(rows,1)
    sid = rows{i,1}; scene_id = rows{i,2}; block_id = rows{i,3}; wwr = rows{i,4}; cx = rows{i,5}; band = rows{i,6}; vec = rows{i,7};
    for ch = 1:nbchan
        out(end+1,:) = {sid, scene_id, block_id, wwr, cx, ch, band, double(vec(ch))}; %#ok<AGROW>
    end
end
T = cell2table(out, 'VariableNames', {'subject_id','scene_id','block_id','WWR','Complexity','chan_idx','band','value'});
writetable(T, fp);
end

function ord = normalize_wwr_local(x)
ord = string(x); ord = strtrim(ord);
for i=1:numel(ord)
    tok = regexp(char(ord(i)), '(\d+)', 'tokens', 'once');
    if ~isempty(tok), ord(i) = string(str2double(tok{1})); end
end
ok = ismember(ord,["15","45","75"]); ord(~ok) = "";
end

function cx = normalize_complexity_local(x)
cx = string(x); cx = strtrim(cx); cl = lower(cx); out = repmat("", numel(cx), 1);
out(ismember(cl,["low","0","c0","complexitylow"])) = "ComplexityLow";
out(ismember(cl,["high","1","c1","complexityhigh"])) = "ComplexityHigh";
out(cx=="ComplexityLow") = "ComplexityLow"; out(cx=="ComplexityHigh") = "ComplexityHigh";
isNum = ~isnan(str2double(cx)); out(isNum & str2double(cx)==0) = "ComplexityLow"; out(isNum & str2double(cx)==1) = "ComplexityHigh";
cx = out;
end

function blk = infer_block_from_scene_id(scene_id)
blk = nan(numel(scene_id),1);
try
    s = double(scene_id);
    blk(s>=1 & s<=6) = 1;
    blk(s>=7 & s<=12) = 2;
catch
end
end

function orderTbl = make_scene_grid_order(meta)
wwrLevels = ["15","45","75"];
cxLevels = ["ComplexityLow","ComplexityHigh"];
rows = {};
meta.WWRn = normalize_wwr_local(meta.WWRn);
meta.CX = normalize_complexity_local(meta.CX);
hasSceneId = ismember('scene_id', meta.Properties.VariableNames);
if ~hasSceneId
    meta.scene_id = nan(height(meta), 1);
end
try
    meta = sortrows(meta, {'CX','WWRn','scene_id'});
catch
end
for r = 1:numel(cxLevels)
    for c = 1:numel(wwrLevels)
        idx = find(string(meta.WWRn)==wwrLevels(c) & string(meta.CX)==cxLevels(r), 1, 'first');
        if ~isempty(idx)
            rows(end+1,:) = {double(meta.scene_id(idx)), string(meta.WWRn(idx)), string(meta.CX(idx))}; %#ok<AGROW>
        else
            rows(end+1,:) = {nan, wwrLevels(c), cxLevels(r)}; %#ok<AGROW>
        end
    end
end
orderTbl = cell2table(rows, 'VariableNames', {'scene_id','WWRn','CX'});
end

function T = add_layout_indices(orderTbl, block_id)
T = orderTbl;
if nargin < 2 || ~isfinite(double(block_id))
    T.block_id = repmat(NaN, height(T), 1);
else
    T.block_id = repmat(double(block_id), height(T), 1);
end
row_idx = nan(height(T),1);
col_idx = nan(height(T),1);
for i = 1:height(T)
    if string(T.CX(i))=="ComplexityLow"
        row_idx(i) = 1;
    elseif string(T.CX(i))=="ComplexityHigh"
        row_idx(i) = 2;
    end
    if string(T.WWRn(i))=="15"
        col_idx(i) = 1;
    elseif string(T.WWRn(i))=="45"
        col_idx(i) = 2;
    elseif string(T.WWRn(i))=="75"
        col_idx(i) = 3;
    end
end
T.row_idx = row_idx;
T.col_idx = col_idx;
T = movevars(T, {'block_id','row_idx','col_idx'}, 'Before', 1);
end

function [comboTopo, comboMeta] = collapse_scene_topos_across_blocks(sceneTopo, meta)
comboTopo = {};
comboMeta = table();
if isempty(sceneTopo) || isempty(meta)
    return;
end
baseMeta = unique(meta(:, {'WWRn','CX'}), 'rows');
baseMeta = baseMeta(strlength(baseMeta.WWRn)>0 & strlength(baseMeta.CX)>0, :);
if isempty(baseMeta)
    return;
end
comboTopo = cell(height(baseMeta), 1);
for i = 1:height(baseMeta)
    idx = find(meta.WWRn==baseMeta.WWRn(i) & meta.CX==baseMeta.CX(i));
    vecs = [];
    for j = idx(:)'
        if j <= numel(sceneTopo) && ~isempty(sceneTopo{j})
            vecs = [vecs, double(sceneTopo{j}(:))]; %#ok<AGROW>
        end
    end
    if ~isempty(vecs)
        comboTopo{i} = mean(vecs, 2, 'omitnan');
    else
        comboTopo{i} = [];
    end
end
comboMeta = baseMeta;
end

function s = short_cx(cx)
s = string(cx);
s(s=="ComplexityLow") = "Low";
s(s=="ComplexityHigh") = "High";
end

function subjRoot = resolve_subjects_root_from_batch(fp_batch)
subjRoot = '';
if nargin < 1 || strlength(string(fp_batch))==0
    return;
end
cur = char(string(fp_batch));
for k = 1:6
    cand = fullfile(cur, 'subjects');
    if exist(cand,'dir')
        subjRoot = cand;
        return;
    end
    [parent, leaf] = fileparts(cur);
    if isempty(parent) || strcmp(parent, cur)
        break;
    end
    % Canonical staged layout: subjects live beside batch under runs/current/
    if strcmpi(leaf, 'batch')
        cand2 = fullfile(parent, 'subjects');
        if exist(cand2,'dir')
            subjRoot = cand2;
            return;
        end
    end
    cur = parent;
end
end

function sid = canonical_subject_id_local(x)
sid = strtrim(string(x));
sid = replace(sid, "\\", "/");
for i=1:numel(sid)
    parts = split(sid(i), "/");
    sid(i) = parts(end);
end
sid = regexprep(sid, '\\.set$', '', 'ignorecase');
sid = strtrim(sid);
end
