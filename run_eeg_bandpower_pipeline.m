function run_eeg_bandpower_pipeline(input_path, config_path)
% RUN_EEG_BANDPOWER_PIPELINE
% VR 场景观看实验 EEG 频段功率分析流水线
%
% 状态机分段规则：
%   1→2: adapt（VR适应）| 2→3: transition（过渡）| 3→4: eyes_closed（闭眼基线）
%   4→7: eyes_open（睁眼基线）| 7→8: view（观看场景）
%   8→9: questionnaire_small（小问卷）| 8→5/6: questionnaire_big（大问卷）
%   9→7/8: gray_to_next | 9→5: gray_to_rest | 9→6: gray_to_end
%   5→7: rest（组间休息）
%
% 输出文件：
%   CSV: bandpower_roi, bandpower_summary, bandpower_tests, scene_level, pairs_check, qc
%   图表: roi_front/par/occ, paired_alpha_recovery, occ_psd, topoplot_theta/alpha/beta
%         paired_scatter, block_comparison, qc_distributions, beta_split, 
%         three_stage_chain, scene_sequence


%% ===== 0) 输入与配置 =====
% 支持：
%   - 不传参：弹窗选择 .set
%   - 传文件路径：处理单个 .set
%   - 传文件夹路径：批量处理该目录下全部 .set
%   - 可选传入 config.json 路径覆盖参数

if nargin < 1 || isempty(input_path)
    input_path = '';
end
if nargin < 2 || isempty(config_path)
    % 默认同目录下的 config.json
    config_path = fullfile(fileparts(mfilename('fullpath')), 'config.json');
end

cfg = load_cfg(config_path);

% 解析输入路径
files = {};
if isempty(input_path)
    choice = questdlg('Select input type', 'EEG Pipeline', 'File', 'Folder', 'Cancel', 'File');
    if strcmp(choice,'File')
        [fn, fp] = uigetfile('*.set', 'Select .set file');
        if isequal(fn,0)
            error('No file selected');
        end
        files = {fullfile(fp, fn)};
    elseif strcmp(choice,'Folder')
        fp = uigetdir(pwd, 'Select folder containing .set files');
        if isequal(fp,0)
            error('No folder selected');
        end
        d = dir(fullfile(fp, '*.set'));
        if isempty(d)
            error('No .set files found in folder: %s', fp);
        end
        files = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
    else
        error('No input selected');
    end
elseif isfolder(input_path)
    d = dir(fullfile(input_path, '*.set'));
    if isempty(d)
        error('No .set files found in folder: %s', input_path);
    end
    files = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
else
    if ~endsWith(lower(input_path), '.set')
        error('Input file must be a .set file');
    end
    files = {input_path};
end

% 可选日志输出
if isfield(cfg, 'log_file') && ~isempty(cfg.log_file)
    try
        diary(cfg.log_file);
    catch
    end
end

%% ===== 0.5) 批量处理 =====
for fi = 1:numel(files)
    this_file = files{fi};
    [fp, name, ext] = fileparts(this_file);
    fn = [name ext];
    fprintf('
=== Processing %s (%d/%d) ===
', fn, fi, numel(files));

    EEG = pop_loadset('filename', fn, 'filepath', fp);
    EEG = eeg_checkset(EEG);

    fs = EEG.srate;

%% ===== 1) 频段、ROI、Welch 参数 =====
bands.theta = [4 7];
bands.alpha = [8 12];
bands.beta  = [13 30];
bands.low_beta  = [13 20];  % low-beta（较少受肌电影响）
bands.high_beta = [20 30];  % high-beta（容易受肌电影响）
totalBand40 = [1 40];  % 原始分母（1-40Hz）
totalBand30 = [1 30];  % 推荐分母（1-30Hz，避免高频噪声影响）

%% ===== 1.5) 可配置参数 (QC 阈值 / 配对模式) =====
if ~exist('cfg', 'var') || isempty(cfg)
    cfg = struct();
end
% QC 阈值：灰屏时长合理范围（秒）
if ~isfield(cfg,'gray_dur_min'); cfg.gray_dur_min = 3; end     % 太短可能是 marker 抖动
if ~isfield(cfg,'gray_dur_max'); cfg.gray_dur_max = 15; end    % 太长可能是异常
% QC 阈值：问卷时长合理范围（秒）
if ~isfield(cfg,'quest_dur_min'); cfg.quest_dur_min = 5; end    % 太短可能是漏填/误触发
if ~isfield(cfg,'quest_dur_max'); cfg.quest_dur_max = 120; end  % 太长可能是异常
% 配对模式：'strict' = 多 gray 时排除该对; 'lenient' = 多 gray 取第一个
if ~isfield(cfg,'pairing_mode'); cfg.pairing_mode = 'strict'; end  % 推荐 strict，避免异常数据污染主结果
% 输出详细日志
if ~isfield(cfg,'verbose'); cfg.verbose = true; end

labels = upper(string({EEG.chanlocs.labels}));
roi.front = find(ismember(labels, ["F3","F4"]));
roi.par   = find(ismember(labels, ["P3","PZ","P4"]));
roi.occ   = find(ismember(labels, ["O1","OZ","O2"]));
assert(~isempty(roi.occ), 'Occipital channels not found. Check channel labels.');

% Welch（2秒窗）
wlen  = round(fs*2);
nover = round(wlen/2);
nfft  = 2^nextpow2(wlen);

%% ===== 2) 状态机分段逻辑 =====
% 核心规则（完全按实验流程）：
% - 1→2 = adapt
% - 3→4 = eyes_closed  
% - 4→7 = eyes_open
% - 7→8 = view（每次进入新 scene）
% - 8→9 = questionnaire_small（正常循环）
% - 8→5 或 8→6 = questionnaire_big（block 结束，跳过灰屏）
% - 9→7 = gray（继续下一个 scene）
% - 9→5 或 9→6 = gray（block 结束前的灰屏）
% - 5→7 = rest（进入下一个 block）

% ---- 收集 marker 事件 (type + latency) ----
types_raw = {EEG.event.type};
lat_raw   = [EEG.event.latency];

mk = nan(size(types_raw));
for i = 1:numel(types_raw)
    if isnumeric(types_raw{i})
        mk(i) = types_raw{i};
    else
        v = str2double(string(types_raw{i}));
        if ~isnan(v), mk(i) = v; end
    end
end

% 只保留 1-9 的有效 marker
keep = ismember(mk, 1:9);
mk   = mk(keep);
lat  = lat_raw(keep);

% 按 latency 排序
[lat, ord] = sort(lat);
mk = mk(ord);

% ---- 状态机分段 ----
nMk = numel(mk);
nSeg = nMk - 1;
segs = struct('seg_idx',[],'m0',[],'m1',[],'s0',[],'s1',[],...
              'start_s',[],'end_s',[],'dur_s',[],'cond',"",...
              'block_id',[],'cycle_in_block',[],'scene_id',[]);

block_id = 1;
scene_in_block = 0;  % block 内的 scene 计数
scene_id_global = 0; % 全局 scene_id
seg_count = 0;

for i = 1:nSeg
    m0 = mk(i);
    m1 = mk(i+1);

    s0 = round(lat(i));
    s1 = round(lat(i+1));

    start_s = (s0-1)/fs;
    end_s   = (s1-1)/fs;
    dur_s   = end_s - start_s;

    cond = "unknown";
    gray_subtype = "";  % gray 子类型
    is_valid_transition = true;  % 转移合法性

    % === 允许的转移表 ===
    % 1→2, 3→4, 4→7, 7→8, 8→9, 8→5, 8→6, 9→7, 9→8, 9→5, 9→6, 5→7
    
    % === 状态机规则 ===
    if m0==1 && m1==2
        cond = "adapt";
        
    elseif m0==2 && m1==3
        % adapt 结束后到基线开始前的过渡段（VR/实验员操作）
        cond = "transition";
    elseif m0==3 && m1==4
        cond = "eyes_closed";
        
    elseif m0==4 && m1==7
        cond = "eyes_open";
        
    elseif m0==7 && m1==8
        cond = "view";
        scene_in_block = scene_in_block + 1;
        scene_id_global = scene_id_global + 1;
        
    elseif m0==8 && m1==9
        cond = "questionnaire_small";
        
    elseif m0==8 && (m1==5 || m1==6)
        cond = "questionnaire_big";
        
    elseif m0==9 && ismember(m1, [7, 8, 5, 6])
        % 9→X 灰屏，只允许 7/8/5/6 作为后继
        cond = "gray";
        % 细分子类型
        if m1==7 || m1==8
            gray_subtype = "gray_to_next";  % 灰屏后进入下一轮
        elseif m1==5
            gray_subtype = "gray_to_rest";  % 灰屏后进入休息
        elseif m1==6
            gray_subtype = "gray_to_end";   % 灰屏后结束
        end
        
    elseif m0==5 && m1==7
        cond = "rest";
        block_id = block_id + 1;
        scene_in_block = 0;
        
    else
        % 非法转移
        cond = "INVALID";
        is_valid_transition = false;
    end

    % 边界检查
    if s0 < 1 || s1 > EEG.pnts || s1 <= s0
        continue;
    end

    % 段长度检查（太短的跳过，但灰屏/休息/适应/过渡段保留）
    minLen = 1.0; % 秒
    if dur_s < minLen
        if ~ismember(cond, ["gray", "rest", "adapt", "transition", "INVALID"])
            continue;
        end
    end

    seg_count = seg_count + 1;
    segs(seg_count).seg_idx = seg_count;
    segs(seg_count).m0 = m0;
    segs(seg_count).m1 = m1;
    segs(seg_count).s0 = s0;
    segs(seg_count).s1 = s1;
    segs(seg_count).start_s = start_s;
    segs(seg_count).end_s = end_s;
    segs(seg_count).dur_s = dur_s;
    segs(seg_count).cond = cond;
    segs(seg_count).gray_subtype = gray_subtype;
    segs(seg_count).is_valid = is_valid_transition;
    segs(seg_count).block_id = block_id;
    segs(seg_count).cycle_in_block = scene_in_block;
    if cond == "view"
        segs(seg_count).scene_id = scene_id_global;
    else
        segs(seg_count).scene_id = NaN;
    end
end

% 移除未使用的 struct 元素
segs = segs(1:seg_count);

assert(~isempty(segs), 'No valid segments extracted. Check marker logic or data.');

% ========== 转移合法性验证 ==========
conds = string({segs.cond});
invalid_segs = find(conds == "INVALID");
if ~isempty(invalid_segs)
    fprintf('\n[ERROR] %d INVALID transitions detected:\n', numel(invalid_segs));
    for ii = invalid_segs
        fprintf('  seg_idx=%d, time=%.1fs, transition=%d→%d\n', ...
            segs(ii).seg_idx, segs(ii).start_s, segs(ii).m0, segs(ii).m1);
    end
else
    fprintf('\n[OK] All marker transitions are valid.\n');
end

% ========== 输出分段统计 ==========
fprintf('\n=== State-machine segmentation ===\n');
[uC,~,ic] = unique(conds, 'stable');
cnt = accumarray(ic, 1);
for ii = 1:numel(uC)
    fprintf('%s: %d\n', uC(ii), cnt(ii));
end

% ========== 按 Block 守恒检查 ==========
fprintf('\n=== Per-block segment counts ===\n');
block_ids = [segs.block_id];
for bid = unique(block_ids)
    mask = block_ids == bid;
    block_conds = conds(mask);
    nv = sum(block_conds=="view");
    ns = sum(block_conds=="questionnaire_small");
    nb = sum(block_conds=="questionnaire_big");
    ng = sum(block_conds=="gray");
    fprintf('Block %d: view=%d, q_small=%d, q_big=%d, gray=%d\n', bid, nv, ns, nb, ng);
    
    % 检查是否符合预期 (每 block: 6 view, 6 small, 1 big, 6 gray)
    if nv==6 && ns==6 && nb==1 && ng==6
        fprintf('  [OK] Block %d counts match expected.\n', bid);
    else
        fprintf('  [NOTE] Block %d counts differ from ideal (6,6,1,6).\n', bid);
    end
end

% ========== 全局守恒检查（分级判定） ==========
n_view = sum(conds=="view");
n_small = sum(conds=="questionnaire_small");
n_big = sum(conds=="questionnaire_big");
n_gray = sum(conds=="gray");
n_rest = sum(conds=="rest");
fprintf('\n[Global] view=%d, q_small=%d, q_big=%d, gray=%d, rest=%d\n', n_view, n_small, n_big, n_gray, n_rest);

% 分级判定
data_quality = 'GOOD';  % GOOD / WARN / FAIL

% Hard fail: view≠12 或 q_big≠2（核心结构错误）
if n_view ~= 12 || n_big ~= 2
    data_quality = 'FAIL';
    fprintf('[HARD FAIL] Critical structure error: view=%d (expected 12), q_big=%d (expected 2)\n', n_view, n_big);
    fprintf('  -> This dataset may be UNUSABLE. Check marker sequence carefully.\n');
% Soft warn: gray≠12 或 q_small≠12（可继续但需注意）
elseif n_gray ~= 12 || n_small ~= 12
    data_quality = 'WARN';
    fprintf('[SOFT WARN] Minor segment count mismatch:\n');
    if n_gray ~= 12
        fprintf('  - gray=%d (expected 12): %d cycles missing gray\n', n_gray, 12-n_gray);
    end
    if n_small ~= 12
        fprintf('  - q_small=%d (expected 12): %d cycles missing questionnaire\n', n_small, 12-n_small);
    end
    fprintf('  -> Analysis will continue. Check pairs_check.csv for affected cycles.\n');
else
    fprintf('[OK] Global counts match expected (12,12,2,12,1)\n');
end

% 保存质量标记到工作区（供后续使用）
segs(1).data_quality = data_quality;

%% ===== 4) 计算每段 ROI 频段功率 + QC指标 =====
Nseg = numel(segs);
% 相对功率（1-30Hz 和 1-40Hz 两个版本）
out_rel30 = nan(Nseg, 15);  % 1-30Hz 版本: 3 ROI x 5 频段 (theta/alpha/beta/low_beta/high_beta)
out_rel40 = nan(Nseg, 15);  % 1-40Hz 版本
out_abs   = nan(Nseg, 15);  % 绝对功率
out_qc    = nan(Nseg, 4);   % QC指标

for i=1:Nseg
    seg = double(EEG.data(:, segs(i).s0:segs(i).s1)); % chan x time
    [P,F] = pwelch(seg', wlen, nover, nfft, fs); % freq x chan

    % ROI 平均PSD后积分（同时计算两个版本的相对功率）
    [out_rel30(i,:), out_rel40(i,:), out_abs(i,:)] = compute_roi_power_v2(P, F, roi, bands, totalBand30, totalBand40);
    
    % QC 指标计算
    out_qc(i,:) = compute_qc_metrics(seg, P, F, segs(i), EEG.pnts, fs);
end

%% ===== 5) 生成表 & 导出CSV =====
T = struct2table(segs);

% === 相对功率 (1-30Hz 分母，推荐使用) ===
% Front
T.F_theta = out_rel30(:,1); T.F_alpha = out_rel30(:,2); T.F_beta = out_rel30(:,3);
T.F_low_beta = out_rel30(:,4); T.F_high_beta = out_rel30(:,5);
% Parietal
T.P_theta = out_rel30(:,6); T.P_alpha = out_rel30(:,7); T.P_beta = out_rel30(:,8);
T.P_low_beta = out_rel30(:,9); T.P_high_beta = out_rel30(:,10);
% Occipital
T.O_theta = out_rel30(:,11); T.O_alpha = out_rel30(:,12); T.O_beta = out_rel30(:,13);
T.O_low_beta = out_rel30(:,14); T.O_high_beta = out_rel30(:,15);

% === 相对功率 (1-40Hz 分母，作为对照) ===
T.F_theta_40 = out_rel40(:,1); T.F_alpha_40 = out_rel40(:,2); T.F_beta_40 = out_rel40(:,3);
T.P_theta_40 = out_rel40(:,6); T.P_alpha_40 = out_rel40(:,7); T.P_beta_40 = out_rel40(:,8);
T.O_theta_40 = out_rel40(:,11); T.O_alpha_40 = out_rel40(:,12); T.O_beta_40 = out_rel40(:,13);

% === log10 绝对功率 ===
T.F_theta_logabs = log10(out_abs(:,1) + eps);
T.F_alpha_logabs = log10(out_abs(:,2) + eps);
T.F_beta_logabs  = log10(out_abs(:,3) + eps);
T.P_theta_logabs = log10(out_abs(:,6) + eps);
T.P_alpha_logabs = log10(out_abs(:,7) + eps);
T.P_beta_logabs  = log10(out_abs(:,8) + eps);
T.O_theta_logabs = log10(out_abs(:,11) + eps);
T.O_alpha_logabs = log10(out_abs(:,12) + eps);
T.O_beta_logabs  = log10(out_abs(:,13) + eps);

% === QC 指标 ===
T.qc_rms_mean   = out_qc(:,1);  % 平均 RMS 振幅
T.qc_hf_ratio   = out_qc(:,2);  % 20-40Hz 高频占比（肌电风险）
T.qc_seg_len_s  = out_qc(:,3);  % 段长度
T.qc_near_boundary = out_qc(:,4); % 是否靠近数据边界

% === 经典比值指标 ===
% TAR: theta/alpha, TBR: theta/beta, BA: beta/alpha
T.F_TAR = T.F_theta ./ T.F_alpha;
T.F_TBR = T.F_theta ./ T.F_beta;
T.F_BA  = T.F_beta  ./ T.F_alpha;

T.P_TAR = T.P_theta ./ T.P_alpha;
T.P_TBR = T.P_theta ./ T.P_beta;
T.P_BA  = T.P_beta  ./ T.P_alpha;

T.O_TAR = T.O_theta ./ T.O_alpha;
T.O_TBR = T.O_theta ./ T.O_beta;
T.O_BA  = T.O_beta  ./ T.O_alpha;

% 条件计数输出
disp('=== Segment counts by condition ===');
conds = string(T.cond);
[uC,~,ic] = unique(conds, 'stable');
cnt = accumarray(ic, 1);
meanDur = accumarray(ic, T.dur_s, [], @mean);
for ii = 1:numel(uC)
    fprintf('%s : %d (mean dur %.1fs)\n', uC(ii), cnt(ii), meanDur(ii));
end

% ===== 添加 pair_id 和 pair_gap_s 审计字段 =====
% 为每个 view 和对应的 gray 分配相同的 pair_id，便于回溯
T.pair_id = nan(height(T), 1);     % 配对编号（view 和 gray 共享）
T.pair_gap_s = nan(height(T), 1);  % 配对间隔（仅 gray 有值）

view_indices = find(conds == "view");
gray_indices = find(conds == "gray");
pair_counter = 0;

for i = 1:numel(view_indices)
    vi = view_indices(i);
    view_end = T.end_s(vi);
    
    if i < numel(view_indices)
        next_view_start = T.start_s(view_indices(i+1));
    else
        next_view_start = Inf;
    end
    
    grays_in_window = gray_indices(T.start_s(gray_indices) > view_end & ...
                                    T.start_s(gray_indices) < next_view_start);
    
    if ~isempty(grays_in_window)
        gi = grays_in_window(1);  % 取第一个 gray（与配对逻辑一致）
        pair_counter = pair_counter + 1;
        T.pair_id(vi) = pair_counter;
        T.pair_id(gi) = pair_counter;
        T.pair_gap_s(gi) = T.start_s(gi) - T.end_s(vi);
    end
end

fprintf('Assigned pair_id to %d view-gray pairs.\n', pair_counter);

% CSV 输出文件名（保证不空）
[base,~,~] = fileparts(fn);
csvfile = fullfile(fp, sprintf('%s_bandpower_roi.csv', base));
writetable(T, csvfile);
fprintf('Saved CSV: %s\n', csvfile);

%% ===== 6) 你要的两行均值 =====
m_closed = mean(T.O_alpha(T.cond=="eyes_closed"), 'omitnan');
m_open   = mean(T.O_alpha(T.cond=="eyes_open"),   'omitnan');
m_view   = mean(T.O_alpha(T.cond=="view"),        'omitnan');
m_gray   = mean(T.O_alpha(T.cond=="gray"),        'omitnan');

fprintf('Occipital relative alpha mean: eyes_closed=%.4f, eyes_open=%.4f\n', m_closed, m_open);
fprintf('Occipital relative alpha mean: view=%.4f, gray=%.4f\n', m_view, m_gray);

%% ===== 7) 可视化图：ROI 条件均值±SEM（θ/α/β） =====
plot_roi_bars(T, fp, base);

%% ===== 8) 可视化图：Occ α view-gray 配对差值（按时间顺序配对） =====
plot_paired_gray_minus_view(T, fp, base);

%% ===== 9) 可视化图：Occ PSD 曲线（各条件均值） =====
plot_occ_psd(EEG, T, roi, fs, wlen, nover, nfft, fp, base);

%% ===== 10) view-gray 地形图（theta/alpha/beta） =====
plot_topoplot_view_minus_gray(EEG, T, bands, totalBand30, fs, wlen, nover, nfft, fp, base);

%% ===== 11) 导出论文级汇总表 =====
export_summary_tables(T, fp, base);

%% ===== 12) 导出 scene_level.csv 用于回归分析 =====
export_scene_level(T, fp, base);

%% ===== 13) 导出配对检查清单（用于审核配对是否正确） =====
export_pairs_check(T, fp, base, cfg);

%% ===== 14) 导出 QC 质量指标表 =====
export_qc_table(T, fp, base);

%% ===== 15) 配对散点+连线图（view vs gray） =====
plot_paired_scatter(T, fp, base, cfg);

%% ===== 16) Block1 vs Block2 稳定性图 =====
plot_block_comparison(T, fp, base);

%% ===== 17) QC 分布图 =====
plot_qc_distributions(T, fp, base, cfg);

%% ===== 18) Low-beta vs High-beta 对比图 =====
plot_beta_split(T, fp, base);

%% ===== 19) 三阶段时间链图（view → quest → gray） =====
plot_three_stage_chain(T, fp, base);

%% ===== 20) Scene 序列曲线 =====
plot_scene_sequence(T, fp, base);

% NOTE: Topoplot with 8 channels is illustrative only (limited spatial resolution).
fprintf('\nNote: Topoplot with %d channels is illustrative only.\n', EEG.nbchan);

disp('Done.');
% end single-file processing
end

end

%% ================== helper functions ==================

function [types, lat_samp, tsec] = get_events_sorted(EEG, fs)
ev = EEG.event;
n = numel(ev);
types = strings(n,1);
lat_samp = zeros(n,1);
tsec = zeros(n,1);

for i=1:n
    typ = ev(i).type;
    if iscell(typ), typ = typ{1}; end
    types(i) = string(typ);
    lat_samp(i) = double(ev(i).latency);
    tsec(i) = (lat_samp(i)-1)/fs;
end

[lat_samp, idx] = sort(lat_samp);
types = types(idx);
tsec  = tsec(idx);
end

function [out_rel30, out_rel40, out_abs] = compute_roi_power_v2(P, F, roi, bands, totalBand30, totalBand40)
% 同时返回两个版本的相对功率（1-30Hz 和 1-40Hz 分母）和绝对功率
% P: freq x chan
% 每个 ROI 输出 5 个频段: theta, alpha, beta, low_beta, high_beta
roiList = {roi.front, roi.par, roi.occ};
out_rel30 = nan(1,15);  % 3 ROI x 5 频段
out_rel40 = nan(1,15);
out_abs   = nan(1,15);
oi = 0;

for r=1:3
    chidx = roiList{r};
    Pm = mean(P(:, chidx), 2, 'omitnan');

    % theta
    oi = oi + 1;
    [out_rel30(oi), out_rel40(oi), out_abs(oi)] = band_power_v2(Pm, F, bands.theta, totalBand30, totalBand40);
    % alpha
    oi = oi + 1;
    [out_rel30(oi), out_rel40(oi), out_abs(oi)] = band_power_v2(Pm, F, bands.alpha, totalBand30, totalBand40);
    % beta (full)
    oi = oi + 1;
    [out_rel30(oi), out_rel40(oi), out_abs(oi)] = band_power_v2(Pm, F, bands.beta, totalBand30, totalBand40);
    % low_beta
    oi = oi + 1;
    [out_rel30(oi), out_rel40(oi), out_abs(oi)] = band_power_v2(Pm, F, bands.low_beta, totalBand30, totalBand40);
    % high_beta
    oi = oi + 1;
    [out_rel30(oi), out_rel40(oi), out_abs(oi)] = band_power_v2(Pm, F, bands.high_beta, totalBand30, totalBand40);
end
end

function [rel30, rel40, absPwr] = band_power_v2(Px, F, band, totalBand30, totalBand40)
% 计算相对功率（两个分母版本）和绝对功率
idx_band = F>=band(1) & F<=band(2);
idx_total30 = F>=totalBand30(1) & F<=totalBand30(2);
idx_total40 = F>=totalBand40(1) & F<=totalBand40(2);

absPwr = trapz(F(idx_band), Px(idx_band));  % 绝对功率 (µV²)
totalPwr30 = trapz(F(idx_total30), Px(idx_total30));
totalPwr40 = trapz(F(idx_total40), Px(idx_total40));

rel30 = absPwr / totalPwr30;  % 相对功率 (1-30Hz 分母)
rel40 = absPwr / totalPwr40;  % 相对功率 (1-40Hz 分母)
end

function [rel, absPwr] = band_power(Px, F, band, totalBand)
% 简化版：计算单一分母的相对功率（供 topoplot 等使用）
idx_band = F>=band(1) & F<=band(2);
idx_total = F>=totalBand(1) & F<=totalBand(2);

absPwr = trapz(F(idx_band), Px(idx_band));
totalPwr = trapz(F(idx_total), Px(idx_total));

rel = absPwr / totalPwr;
end

function qc = compute_qc_metrics(seg, P, F, segInfo, totalPnts, fs)
% 计算 QC 指标
% seg: chan x time
% 返回: [rms_mean, hf_ratio, seg_len_s, near_boundary]

% 1) 平均 RMS 振幅（所有通道）
rms_vals = rms(seg, 2);  % 每通道的 RMS
rms_mean = mean(rms_vals);

% 2) 20-40Hz 高频功率占比（β/肌电风险指标）
% 用所有通道平均 PSD
Pm = mean(P, 2, 'omitnan');
idx_hf = F>=20 & F<=40;
idx_total = F>=1 & F<=40;
hf_pwr = trapz(F(idx_hf), Pm(idx_hf));
total_pwr = trapz(F(idx_total), Pm(idx_total));
hf_ratio = hf_pwr / total_pwr;

% 3) 段长度（秒）
seg_len_s = segInfo.dur_s;

% 4) 是否靠近数据边界（前后 2 秒内）
boundary_margin = 2 * fs;  % 2 秒
near_start = segInfo.s0 < boundary_margin;
near_end   = segInfo.s1 > (totalPnts - boundary_margin);
near_boundary = double(near_start || near_end);

qc = [rms_mean, hf_ratio, seg_len_s, near_boundary];
end

function plot_roi_bars(T, fp, base)
conds = unique(string(T.cond),'stable');
roiNames = ["Front(F3/F4)","Par(P3/Pz/P4)","Occ(O1/Oz/O2)"];
roiShort = ["front","par","occ"];
bandNames = ["theta","alpha","beta"];

% 组织数据：每个cond的均值和SEM
vals = struct();
vals.front = [T.F_theta, T.F_alpha, T.F_beta];
vals.par   = [T.P_theta, T.P_alpha, T.P_beta];
vals.occ   = [T.O_theta, T.O_alpha, T.O_beta];

S = struct();
fields = ["front","par","occ"];

for fi = 1:numel(fields)
    f = fields(fi);
    M = nan(numel(conds),3);
    E = nan(numel(conds),3);
    for i=1:numel(conds)
        idx = string(T.cond)==conds(i);
        X = vals.(f)(idx,:);
        M(i,:) = mean(X, 1, 'omitnan');
        n = sum(~isnan(X(:,1)));
        E(i,:) = std(X,0,1,'omitnan') ./ sqrt(max(n,1));
    end

    fig = figure('Name', sprintf('ROI %s mean ± SEM', f));
    b = bar(M); hold on;
    % 误差条：逐band画（不手动指定颜色）
    for j=1:3
        x = b(j).XEndPoints;
        errorbar(x, M(:,j), E(:,j), '.', 'LineWidth', 1);
    end
    set(gca,'XTickLabel',conds,'XTickLabelRotation',30);
    legend(bandNames,'Location','best');
    ylabel('Relative power');
    title(char(roiNames(fi)));
    grid on;
    
    % 保存图片
    figFile = fullfile(fp, sprintf('%s_roi_%s.png', base, roiShort(fi)));
    saveas(fig, figFile);
    fprintf('Saved figure: %s\n', figFile);
end
end

function plot_paired_gray_minus_view(T, fp, base)
% 就近窗口匹配：对每个 view，在它结束后、下一个 view 开始前，找第一个 gray

conds = string(T.cond);
view_indices = find(conds == "view");
gray_indices = find(conds == "gray");

pairs = [];
for i = 1:numel(view_indices)
    vi = view_indices(i);
    view_end = T.end_s(vi);
    
    if i < numel(view_indices)
        next_view_start = T.start_s(view_indices(i+1));
    else
        next_view_start = Inf;
    end
    
    grays_in_window = gray_indices(T.start_s(gray_indices) > view_end & ...
                                    T.start_s(gray_indices) < next_view_start);
    
    if ~isempty(grays_in_window)
        pairs(end+1,:) = [vi, grays_in_window(1)]; %#ok<AGROW>
    end
end

if isempty(pairs)
    warning('No view-gray pairs found.'); 
    return;
end

viewA = T.O_alpha(pairs(:,1));
grayA = T.O_alpha(pairs(:,2));
dA = grayA - viewA;

fig = figure('Name','Per-cycle alpha recovery (gray - view)');
plot(dA,'-o'); grid on;
xlabel('Cycle #'); ylabel('Delta O\_alpha (gray - view)');
title('Occipital alpha recovery per cycle');

% 保存图片
figFile = fullfile(fp, sprintf('%s_paired_alpha_recovery.png', base));
saveas(fig, figFile);
fprintf('Saved figure: %s\n', figFile);

% 同时对 theta/alpha/beta 做 gray vs view 配对检验（三个 ROI）
if size(pairs,1) >= 6
    roiNames = {'F','P','O'};
    bandNames = {'theta','alpha','beta'};
    
    fprintf('\n=== Paired signrank tests: gray vs view (n=%d pairs) ===\n', size(pairs,1));
    
    for ri = 1:3
        roiN = roiNames{ri};
        for bi = 1:3
            bandN = bandNames{bi};
            colName = sprintf('%s_%s', roiN, bandN);
            
            viewVal = T.(colName)(pairs(:,1));
            grayVal = T.(colName)(pairs(:,2));
            
            pVal = signrank(grayVal, viewVal);
            meanDiff = mean(grayVal - viewVal);
            
            fprintf('%s: p=%.4g, mean_diff=%.4f (gray-view)\n', colName, pVal, meanDiff);
        end
    end
end
end

function plot_occ_psd(EEG, T, roi, fs, wlen, nover, nfft, fp, base)
condsToPlot = unique(string(T.cond),'stable');

fig = figure('Name','Occipital PSD (mean over ROI)');
hold on;

for c = condsToPlot'
    idxc = find(string(T.cond)==c);
    if isempty(idxc), continue; end

    Pall = [];
    for kk = idxc'
        seg = double(EEG.data(:, T.s0(kk):T.s1(kk)));
        [P,F] = pwelch(seg', wlen, nover, nfft, fs); % freq x chan
        if isempty(Pall), Pall = zeros(size(P,1),1); end
        Pall = Pall + mean(P(:, roi.occ), 2, 'omitnan');
    end
    Pall = Pall / numel(idxc);

    plot(F, 10*log10(Pall), 'DisplayName', char(c));
end

xlim([1 40]); xlabel('Hz'); ylabel('Power (dB)');
legend('show'); grid on;
title('Occipital PSD across conditions');

% 保存图片
figFile = fullfile(fp, sprintf('%s_occ_psd.png', base));
saveas(fig, figFile);
fprintf('Saved figure: %s\n', figFile);
end

function plot_topoplot_view_minus_gray(EEG, T, bands, totalBand, fs, wlen, nover, nfft, fp, base)
% 绘制 theta/alpha/beta 三个频段的 view-gray 差异地形图

% 检查 chanlocs
if ~isfield(EEG,'chanlocs') || numel(EEG.chanlocs)~=EEG.nbchan
    warning('No valid chanlocs found. Skipping topoplot.');
    return;
end

idx_view = find(string(T.cond)=="view");
idx_gray = find(string(T.cond)=="gray");
if isempty(idx_view) || isempty(idx_gray)
    warning('No segments for view/gray; skip topoplot.'); 
    return;
end

% 频段列表
bandNames = {'theta', 'alpha', 'beta'};
bandRanges = {bands.theta, bands.alpha, bands.beta};

for bi = 1:3
    bandName = bandNames{bi};
    bandRange = bandRanges{bi};
    
    rB_view = zeros(EEG.nbchan,1); cntv = 0;
    rB_gray = zeros(EEG.nbchan,1); cntg = 0;
    
    % 计算 view 段均值
    for kk = idx_view'
        seg = double(EEG.data(:, T.s0(kk):T.s1(kk)));
        [P,F] = pwelch(seg', round(fs*2), round(fs), 2^nextpow2(round(fs*2)), fs);
        for ci = 1:EEG.nbchan
            [rel, ~] = band_power(P(:,ci), F, bandRange, totalBand);
            rB_view(ci) = rB_view(ci) + rel;
        end
        cntv = cntv + 1;
    end
    rB_view = rB_view / cntv;
    
    % 计算 gray 段均值
    for kk = idx_gray'
        seg = double(EEG.data(:, T.s0(kk):T.s1(kk)));
        [P,F] = pwelch(seg', round(fs*2), round(fs), 2^nextpow2(round(fs*2)), fs);
        for ci = 1:EEG.nbchan
            [rel, ~] = band_power(P(:,ci), F, bandRange, totalBand);
            rB_gray(ci) = rB_gray(ci) + rel;
        end
        cntg = cntg + 1;
    end
    rB_gray = rB_gray / cntg;
    
    % 差异
    dB = rB_view - rB_gray;
    
    % 绘图
    fig = figure('Name', sprintf('%s topoplot (view - gray)', upper(bandName)));
    topoplot(dB, EEG.chanlocs, 'electrodes', 'labels');
    title(sprintf('Relative %s difference: view - gray', bandName));
    colorbar;
    
    % 保存图片
    figFile = fullfile(fp, sprintf('%s_topoplot_%s.png', base, bandName));
    saveas(fig, figFile);
    fprintf('Saved figure: %s\n', figFile);
end
end

function export_summary_tables(T, fp, base)
% 输出一个 summary.csv：每个cond在各ROI各频段的 N/mean/SEM/median

roiList = ["F","P","O"];
bandList = ["theta","alpha","beta"];
conds = unique(string(T.cond),'stable');

rows = {};
for r = 1:numel(roiList)
    for b = 1:numel(bandList)
        col = sprintf('%s_%s', roiList(r), bandList(b)); % 例如 O_alpha
        col = char(col);

        for c = 1:numel(conds)
            idx = string(T.cond)==conds(c);
            x = T.(col)(idx);
            x = x(~isnan(x));

            N = numel(x);
            if N==0
                mu = NaN; sem = NaN; med = NaN;
            else
                mu = mean(x);
                sem = std(x)/sqrt(N);
                med = median(x);
            end

            rows(end+1,:) = {roiList(r), bandList(b), conds(c), N, mu, sem, med}; %#ok<AGROW>
        end
    end
end

S = cell2table(rows, 'VariableNames', ...
    {'ROI','Band','Cond','N','Mean','SEM','Median'});

% 加入两条你关心的对比 p 值（如果条件存在）
p_view_gray = NaN;
p_close_open = NaN;

if any(conds=="view") && any(conds=="gray")
    % 就近窗口匹配
    view_indices = find(conds == "view");
    gray_indices = find(conds == "gray");
    pairs = [];
    for i = 1:numel(view_indices)
        vi = view_indices(i);
        view_end = T.end_s(vi);
        if i < numel(view_indices)
            next_view_start = T.start_s(view_indices(i+1));
        else
            next_view_start = Inf;
        end
        grays_in_window = gray_indices(T.start_s(gray_indices) > view_end & ...
                                        T.start_s(gray_indices) < next_view_start);
        if ~isempty(grays_in_window)
            pairs(end+1,:) = [vi, grays_in_window(1)]; %#ok<AGROW>
        end
    end
    if ~isempty(pairs) && size(pairs,1) >= 6
        viewA = T.O_alpha(pairs(:,1));
        grayA = T.O_alpha(pairs(:,2));
        p_view_gray = signrank(grayA, viewA);
    end
end

if any(conds=="eyes_closed") && any(conds=="eyes_open")
    xc = T.O_alpha(string(T.cond)=="eyes_closed");
    xo = T.O_alpha(string(T.cond)=="eyes_open");
    if numel(xc)>=1 && numel(xo)>=1
        % 只有1段时p值没意义，仍保留NaN；这里只记录差值
        p_close_open = NaN;
    end
end

% 把这些写进一个小表附在后面（或单独保存）
meta = table( ...
    ["Occ_alpha_gray_vs_view_p"; "Occ_alpha_eyesClosed_vs_eyesOpen_p"], ...
    [p_view_gray; p_close_open], ...
    'VariableNames', {'Test','Value'});

% 保存
summaryFile = fullfile(fp, sprintf('%s_bandpower_summary.csv', base));
metaFile    = fullfile(fp, sprintf('%s_bandpower_tests.csv', base));

writetable(S, summaryFile);
writetable(meta, metaFile);

fprintf('Saved summary: %s\n', summaryFile);
fprintf('Saved tests:   %s\n', metaFile);
end

function export_scene_level(T, fp, base)
% 导出 scene_level.csv：仅 view 段，按 scene_id 排序
% 包含 scene_id, 各频段功率，可直接用于与复杂度指标回归

% 筛选 view 段
view_mask = string(T.cond) == "view";
Tv = T(view_mask, :);

if isempty(Tv)
    warning('No view segments found. Skipping scene_level export.');
    return;
end

% 按 scene_id 排序
[~, ord] = sort(Tv.scene_id);
Tv = Tv(ord, :);

% 选择回归分析需要的列
sceneTable = table( ...
    Tv.scene_id, ...
    Tv.dur_s, ...
    Tv.F_theta, Tv.F_alpha, Tv.F_beta, ...
    Tv.P_theta, Tv.P_alpha, Tv.P_beta, ...
    Tv.O_theta, Tv.O_alpha, Tv.O_beta, ...
    'VariableNames', { ...
        'scene_id', 'duration_s', ...
        'F_theta', 'F_alpha', 'F_beta', ...
        'P_theta', 'P_alpha', 'P_beta', ...
        'O_theta', 'O_alpha', 'O_beta' ...
    });

% 保存
sceneFile = fullfile(fp, sprintf('%s_scene_level.csv', base));
writetable(sceneTable, sceneFile);
fprintf('Saved scene_level: %s\n', sceneFile);

% 打印预览
disp('=== Scene-level data (for regression) ===');
disp(sceneTable);
end

function export_pairs_check(T, fp, base, cfg)
% 导出配对检查清单：就近窗口匹配
% 配对关系说明：
%   view (7→8) → questionnaire_small (8→9) → gray (9→7/8/5/6)
%   gray 表征 post-view recovery（观看后恢复段），其前通常有问卷阶段
%   因此 gap_sec = gray_start - view_end 包含了问卷填写时间
%
% 配对模式（cfg.pairing_mode）：
%   'strict' = 多 gray 时排除该对，不纳入主统计
%   'lenient' = 多 gray 时取第一个（用于补充分析/敏感性分析）

nRows = height(T);
conds = string(T.cond);

view_indices = find(conds == "view");
gray_indices = find(conds == "gray");

pairs = [];
pair_status = {};  % 'normal' / 'no_gray' / 'multi_gray'
anomalies = {};

for i = 1:numel(view_indices)
    vi = view_indices(i);
    view_end = T.end_s(vi);
    
    if i < numel(view_indices)
        next_view_start = T.start_s(view_indices(i+1));
    else
        next_view_start = Inf;
    end
    
    grays_in_window = gray_indices(T.start_s(gray_indices) > view_end & ...
                                    T.start_s(gray_indices) < next_view_start);
    
    if isempty(grays_in_window)
        % 无 gray：异常
        pair_status{end+1} = 'no_gray'; %#ok<AGROW>
        pairs(end+1,:) = [vi, NaN]; %#ok<AGROW>
        anomalies{end+1} = sprintf('view seg_idx=%d (scene=%d): NO GRAY found', ...
            vi, T.scene_id(vi)); %#ok<AGROW>
            
    elseif numel(grays_in_window) > 1
        % 多 gray：根据模式处理
        if strcmp(cfg.pairing_mode, 'strict')
            % strict 模式：排除该对
            pair_status{end+1} = 'multi_gray_excluded'; %#ok<AGROW>
            pairs(end+1,:) = [vi, NaN]; %#ok<AGROW>
            anomalies{end+1} = sprintf('view seg_idx=%d (scene=%d): MULTIPLE GRAYS (%d), EXCLUDED in strict mode', ...
                vi, T.scene_id(vi), numel(grays_in_window)); %#ok<AGROW>
        else
            % lenient 模式：取第一个
            gi = grays_in_window(1);
            pair_status{end+1} = 'multi_gray_first'; %#ok<AGROW>
            pairs(end+1,:) = [vi, gi]; %#ok<AGROW>
            anomalies{end+1} = sprintf('view seg_idx=%d (scene=%d): MULTIPLE GRAYS (%d), using first (lenient)', ...
                vi, T.scene_id(vi), numel(grays_in_window)); %#ok<AGROW>
        end
    else
        % 正常：刚好一个 gray
        pair_status{end+1} = 'normal'; %#ok<AGROW>
        pairs(end+1,:) = [vi, grays_in_window(1)]; %#ok<AGROW>
    end
end

% 统计有效配对数
valid_pairs = ~isnan(pairs(:,2));
n_valid = sum(valid_pairs);
n_total = size(pairs, 1);

fprintf('\n=== Pairs Check (mode=%s) ===\n', cfg.pairing_mode);
fprintf('Total view segments: %d\n', n_total);
fprintf('Valid pairs: %d (%.0f%%)\n', n_valid, 100*n_valid/n_total);

if n_valid == 0
    warning('No valid view-gray pairs found. Skipping pairs_check export.');
    return;
end

% 构建检查表（仅有效配对）
valid_idx = find(valid_pairs);
view_idx = pairs(valid_idx, 1);
gray_idx = pairs(valid_idx, 2);
nPairs = numel(view_idx);

pairTable = table( ...
    (1:nPairs)', ...
    T.scene_id(view_idx), ...
    T.block_id(view_idx), ...
    T.cycle_in_block(view_idx), ...
    string(pair_status(valid_idx))', ...  % 配对状态
    T.m0(view_idx), T.m1(view_idx), ...
    T.m0(gray_idx), T.m1(gray_idx), ...
    T.gray_subtype(gray_idx), ...
    T.start_s(view_idx), ...
    T.end_s(view_idx), ...
    T.dur_s(view_idx), ...
    T.start_s(gray_idx), ...
    T.end_s(gray_idx), ...
    T.dur_s(gray_idx), ...
    T.start_s(gray_idx) - T.end_s(view_idx), ...
    T.O_alpha(view_idx), ...
    T.O_alpha(gray_idx), ...
    T.O_alpha(gray_idx) - T.O_alpha(view_idx), ...
    T.O_theta(view_idx), ...
    T.O_theta(gray_idx), ...
    T.O_beta(view_idx), ...
    T.O_beta(gray_idx), ...
    'VariableNames', { ...
        'pair_id', 'scene_id', 'block_id', 'cycle_in_block', 'pair_status', ...
        'view_m0', 'view_m1', 'gray_m0', 'gray_m1', 'gray_subtype', ...
        'view_start', 'view_end', 'view_dur', ...
        'gray_start', 'gray_end', 'gray_dur', 'gap_sec', ...
        'view_O_alpha', 'gray_O_alpha', 'delta_O_alpha', ...
        'view_O_theta', 'gray_O_theta', 'view_O_beta', 'gray_O_beta' ...
    });

% 保存
pairsFile = fullfile(fp, sprintf('%s_pairs_check.csv', base));
writetable(pairTable, pairsFile);
fprintf('Saved pairs_check: %s\n', pairsFile);

% 打印摘要
fprintf('\n=== Valid Pairs Summary (%d pairs) ===\n', nPairs);
fprintf('Gap includes questionnaire time (view→quest→gray)\n');
fprintf('Mean gap (view_end to gray_start): %.2f sec\n', mean(pairTable.gap_sec));
fprintf('Mean view duration: %.1f sec\n', mean(pairTable.view_dur));
fprintf('Mean gray duration: %.1f sec\n', mean(pairTable.gray_dur));
fprintf('Mean delta_O_alpha (gray-view): %.4f\n', mean(pairTable.delta_O_alpha));

% 时长合理性检查（使用 cfg 阈值）
short_grays = pairTable.gray_dur < cfg.gray_dur_min;
long_grays = pairTable.gray_dur > cfg.gray_dur_max;
if any(short_grays)
    fprintf('\n[NOTE] %d grays have unusually short duration (<%ds)\n', sum(short_grays), cfg.gray_dur_min);
end
if any(long_grays)
    fprintf('[NOTE] %d grays have unusually long duration (>%ds)\n', sum(long_grays), cfg.gray_dur_max);
end

% 报告异常
if ~isempty(anomalies)
    fprintf('\n[WARNING] %d pairing anomalies detected:\n', numel(anomalies));
    for i = 1:numel(anomalies)
        fprintf('  %s\n', anomalies{i});
    end
else
    fprintf('\n[OK] All view-gray pairs matched normally.\n');
end
end

function export_qc_table(T, fp, base)
% 导出 QC 质量指标表
% 包含每段的 RMS、高频占比、边界标记等

qcTable = table( ...
    T.seg_idx, ...
    string(T.cond), ...
    T.scene_id, ...
    T.start_s, ...
    T.end_s, ...
    T.dur_s, ...
    T.qc_rms_mean, ...
    T.qc_hf_ratio, ...
    T.qc_near_boundary, ...
    'VariableNames', { ...
        'seg_idx', 'cond', 'scene_id', 'start_s', 'end_s', 'dur_s', ...
        'rms_mean_uV', 'hf_ratio_20_40Hz', 'near_boundary' ...
    });

% 保存
qcFile = fullfile(fp, sprintf('%s_qc.csv', base));
writetable(qcTable, qcFile);
fprintf('Saved QC: %s\n', qcFile);

% 打印 QC 摘要
fprintf('\n=== QC Summary ===\n');

% 按条件输出平均值
conds = unique(string(T.cond),'stable');
for c = conds'
    idx = string(T.cond)==c;
    fprintf('%s: mean_RMS=%.2f µV, mean_HF_ratio=%.3f\n', ...
        c, mean(T.qc_rms_mean(idx),'omitnan'), mean(T.qc_hf_ratio(idx),'omitnan'));
end

% 检查高频异常（HF_ratio > 0.4 可能有肌电问题）
high_hf = T.qc_hf_ratio > 0.4;
if any(high_hf)
    fprintf('\n[WARNING] %d segments have HF_ratio > 0.4 (possible EMG contamination)\n', sum(high_hf));
    high_hf_segs = T(high_hf, {'seg_idx','cond','qc_hf_ratio'});
    disp(high_hf_segs);
else
    fprintf('\n[OK] No segments with high-frequency anomaly detected.\n');
end

% 检查边界问题
near_bound = T.qc_near_boundary > 0;
if any(near_bound)
    fprintf('[NOTE] %d segments are near data boundaries.\n', sum(near_bound));
end
end

%% ================== 高级可视化函数 ==================

function plot_paired_scatter(T, fp, base, cfg)
% 配对散点+连线图：view vs gray
% 让配对差异一眼可信，能看出是否被离群点"带出来的"

conds = string(T.cond);
view_indices = find(conds == "view");
gray_indices = find(conds == "gray");

% 构建配对
pairs = [];
for i = 1:numel(view_indices)
    vi = view_indices(i);
    view_end = T.end_s(vi);
    if i < numel(view_indices)
        next_view_start = T.start_s(view_indices(i+1));
    else
        next_view_start = Inf;
    end
    grays_in_window = gray_indices(T.start_s(gray_indices) > view_end & ...
                                    T.start_s(gray_indices) < next_view_start);
    if ~isempty(grays_in_window) && numel(grays_in_window) == 1
        pairs(end+1,:) = [vi, grays_in_window(1)]; %#ok<AGROW>
    end
end

if isempty(pairs)
    warning('No valid pairs for scatter plot.');
    return;
end

% 绘制关键 ROI x Band 组合
plots_to_make = {'O_alpha', 'F_theta', 'O_theta', 'F_alpha'};
titles = {'Occipital Alpha', 'Frontal Theta', 'Occipital Theta', 'Frontal Alpha'};

fig = figure('Name', 'Paired Scatter: View vs Gray', 'Position', [100 100 1000 800]);

for pi = 1:numel(plots_to_make)
    colName = plots_to_make{pi};
    if ~ismember(colName, T.Properties.VariableNames)
        continue;
    end
    
    viewVals = T.(colName)(pairs(:,1));
    grayVals = T.(colName)(pairs(:,2));
    
    subplot(2, 2, pi);
    hold on;
    
    % 绘制连线
    for i = 1:numel(viewVals)
        plot([viewVals(i), grayVals(i)], [1, 2], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end
    
    % 绘制散点
    scatter(viewVals, ones(size(viewVals)), 50, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
    scatter(grayVals, 2*ones(size(grayVals)), 50, 'r', 'filled', 'MarkerFaceAlpha', 0.7);
    
    % 绘制均值
    plot([mean(viewVals,'omitnan'), mean(grayVals,'omitnan')], [1, 2], 'k-', 'LineWidth', 2);
    
    set(gca, 'YTick', [1 2], 'YTickLabel', {'View', 'Gray'});
    xlabel(strrep(colName, '_', ' '));
    title(titles{pi});
    
    % 计算统计
    diff_vals = grayVals - viewVals;
    median_diff = median(diff_vals, 'omitnan');
    if numel(viewVals) >= 6
        pVal = signrank(grayVals, viewVals);
    else
        pVal = NaN;
    end
    
    % 标注统计信息
    text(min(xlim), 2.3, sprintf('n=%d, med.diff=%.3f, p=%.3g', numel(viewVals), median_diff, pVal), ...
        'FontSize', 9);
    
    ylim([0.5 2.8]);
    grid on;
end

sgtitle('Paired View-Gray Comparison (per scene)', 'FontWeight', 'bold');

figFile = fullfile(fp, sprintf('%s_paired_scatter.png', base));
saveas(fig, figFile);
fprintf('Saved figure: %s\n', figFile);
end

function plot_block_comparison(T, fp, base)
% Block1 vs Block2 稳定性图
% 检查疲劳/适应导致的漂移

conds = string(T.cond);
view_mask = conds == "view";
Tv = T(view_mask, :);

if isempty(Tv) || numel(unique(Tv.block_id)) < 2
    warning('Not enough blocks for comparison. Skipping.');
    return;
end

% 关键指标
metrics = {'O_alpha', 'O_theta', 'F_alpha', 'F_theta'};
titles = {'Occ Alpha', 'Occ Theta', 'Front Alpha', 'Front Theta'};

fig = figure('Name', 'Block Comparison', 'Position', [100 100 1000 400]);

for mi = 1:numel(metrics)
    colName = metrics{mi};
    if ~ismember(colName, Tv.Properties.VariableNames)
        continue;
    end
    
    b1_vals = Tv.(colName)(Tv.block_id == 1);
    b2_vals = Tv.(colName)(Tv.block_id == 2);
    
    subplot(1, 4, mi);
    hold on;
    
    % 箱线图
    boxplot([b1_vals; b2_vals], [ones(size(b1_vals)); 2*ones(size(b2_vals))], ...
        'Labels', {'Block 1', 'Block 2'});
    
    % 均值点
    plot(1, mean(b1_vals,'omitnan'), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(2, mean(b2_vals,'omitnan'), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    ylabel(strrep(colName, '_', ' '));
    title(titles{mi});
    
    % 配对检验
    if numel(b1_vals) == numel(b2_vals) && numel(b1_vals) >= 3
        pVal = signrank(b1_vals, b2_vals);
        text(1.5, max([b1_vals; b2_vals])*1.05, sprintf('p=%.3g', pVal), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    
    grid on;
end

sgtitle('Block 1 vs Block 2 Stability (view segments)', 'FontWeight', 'bold');

figFile = fullfile(fp, sprintf('%s_block_comparison.png', base));
saveas(fig, figFile);
fprintf('Saved figure: %s\n', figFile);
end

function plot_qc_distributions(T, fp, base, cfg)
% QC 分布图：让 β/高频解释站得住

fig = figure('Name', 'QC Distributions', 'Position', [100 100 1200 400]);

% 1) HF_ratio 按条件分组
subplot(1, 3, 1);
conds = string(T.cond);
conds_to_plot = ["view", "gray", "questionnaire_small"];
hf_data = [];
hf_groups = [];
for ci = 1:numel(conds_to_plot)
    mask = conds == conds_to_plot(ci);
    vals = T.qc_hf_ratio(mask);
    hf_data = [hf_data; vals]; %#ok<AGROW>
    hf_groups = [hf_groups; ci * ones(size(vals))]; %#ok<AGROW>
end
if ~isempty(hf_data)
    boxplot(hf_data, hf_groups, 'Labels', conds_to_plot);
    ylabel('HF Ratio (20-40Hz / 1-40Hz)');
    title('High-Frequency Contamination by Condition');
    yline(0.4, 'r--', 'LineWidth', 1.5);  % 阈值线
    text(0.5, 0.42, 'EMG threshold', 'Color', 'r', 'FontSize', 8);
end
grid on;

% 2) pair_gap_s 直方图
subplot(1, 3, 2);
gap_vals = T.pair_gap_s(~isnan(T.pair_gap_s));
if ~isempty(gap_vals)
    histogram(gap_vals, 15, 'FaceColor', [0.3 0.6 0.9]);
    xlabel('Gap (view\_end to gray\_start), sec');
    ylabel('Count');
    title('Pair Gap Distribution');
    xline(mean(gap_vals), 'r-', 'LineWidth', 2);
    text(mean(gap_vals)*1.1, max(ylim)*0.9, sprintf('mean=%.1fs', mean(gap_vals)), 'Color', 'r');
end
grid on;

% 3) RMS vs HF_ratio 散点
subplot(1, 3, 3);
hold on;
colors = lines(numel(conds_to_plot));
for ci = 1:numel(conds_to_plot)
    mask = conds == conds_to_plot(ci);
    scatter(T.qc_rms_mean(mask), T.qc_hf_ratio(mask), 30, colors(ci,:), 'filled', 'MarkerFaceAlpha', 0.6);
end
xlabel('RMS Mean (µV)');
ylabel('HF Ratio');
title('RMS vs HF Ratio');
legend(conds_to_plot, 'Location', 'best');
yline(0.4, 'r--', 'LineWidth', 1);
grid on;

sgtitle('Quality Control Distributions', 'FontWeight', 'bold');

figFile = fullfile(fp, sprintf('%s_qc_distributions.png', base));
saveas(fig, figFile);
fprintf('Saved figure: %s\n', figFile);
end

function plot_beta_split(T, fp, base)
% Low-beta vs High-beta 对比图
% 主动区分肌电风险段

conds_to_plot = ["view", "gray", "eyes_closed", "eyes_open"];
conds = string(T.cond);

fig = figure('Name', 'Beta Split: Low vs High', 'Position', [100 100 1000 400]);

rois = {'F', 'P', 'O'};
roi_names = {'Frontal', 'Parietal', 'Occipital'};

for ri = 1:numel(rois)
    subplot(1, 3, ri);
    
    low_col = sprintf('%s_low_beta', rois{ri});
    high_col = sprintf('%s_high_beta', rois{ri});
    
    if ~ismember(low_col, T.Properties.VariableNames)
        continue;
    end
    
    % 计算每个条件的均值
    low_means = [];
    high_means = [];
    low_sems = [];
    high_sems = [];
    valid_conds = {};
    
    for ci = 1:numel(conds_to_plot)
        mask = conds == conds_to_plot(ci);
        if sum(mask) == 0, continue; end
        
        low_vals = T.(low_col)(mask);
        high_vals = T.(high_col)(mask);
        
        low_means(end+1) = mean(low_vals, 'omitnan'); %#ok<AGROW>
        high_means(end+1) = mean(high_vals, 'omitnan'); %#ok<AGROW>
        low_sems(end+1) = std(low_vals,'omitnan') / sqrt(sum(~isnan(low_vals))); %#ok<AGROW>
        high_sems(end+1) = std(high_vals,'omitnan') / sqrt(sum(~isnan(high_vals))); %#ok<AGROW>
        valid_conds{end+1} = char(conds_to_plot(ci)); %#ok<AGROW>
    end
    
    if isempty(low_means), continue; end
    
    x = 1:numel(valid_conds);
    bar_data = [low_means; high_means]';
    
    b = bar(x, bar_data, 'grouped');
    hold on;
    
    % 误差条
    for bi = 1:2
        x_err = b(bi).XEndPoints;
        if bi == 1
            errorbar(x_err, low_means, low_sems, 'k.', 'LineWidth', 1);
        else
            errorbar(x_err, high_means, high_sems, 'k.', 'LineWidth', 1);
        end
    end
    
    set(gca, 'XTick', x, 'XTickLabel', valid_conds, 'XTickLabelRotation', 30);
    ylabel('Relative Power');
    title(roi_names{ri});
    legend({'Low β (13-20Hz)', 'High β (20-30Hz)'}, 'Location', 'best');
    grid on;
end

sgtitle('Low-Beta vs High-Beta (EMG sensitivity comparison)', 'FontWeight', 'bold');

figFile = fullfile(fp, sprintf('%s_beta_split.png', base));
saveas(fig, figFile);
fprintf('Saved figure: %s\n', figFile);
end

function plot_three_stage_chain(T, fp, base)
% 三阶段时间链图：view → questionnaire_small → gray
% 看变化发生在"离开场景时"还是"灰屏恢复时"

conds = string(T.cond);

% 只用有完整三阶段的数据
% 对于每个 view，找对应的 questionnaire 和 gray

view_indices = find(conds == "view");
chains = [];  % [view_idx, quest_idx, gray_idx]

for i = 1:numel(view_indices)
    vi = view_indices(i);
    if vi+2 > height(T), continue; end
    
    c1 = conds(vi+1);
    c2 = conds(vi+2);
    
    if contains(c1, "questionnaire") && c2 == "gray"
        chains(end+1,:) = [vi, vi+1, vi+2]; %#ok<AGROW>
    end
end

if size(chains,1) < 3
    warning('Not enough complete chains for three-stage plot. Skipping.');
    return;
end

% 绘制关键指标
metrics = {'O_alpha', 'F_theta'};
titles = {'Occipital Alpha', 'Frontal Theta'};

fig = figure('Name', 'Three-Stage Chain', 'Position', [100 100 800 400]);

for mi = 1:numel(metrics)
    colName = metrics{mi};
    if ~ismember(colName, T.Properties.VariableNames)
        continue;
    end
    
    subplot(1, 2, mi);
    hold on;
    
    % 每条链一条线
    for ci = 1:size(chains,1)
        vals = T.(colName)(chains(ci,:));
        plot(1:3, vals, 'o-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'MarkerSize', 4);
    end
    
    % 均值线
    means = zeros(1,3);
    sems = zeros(1,3);
    for si = 1:3
        vals = T.(colName)(chains(:,si));
        means(si) = mean(vals, 'omitnan');
        sems(si) = std(vals,'omitnan') / sqrt(numel(vals));
    end
    
    errorbar(1:3, means, sems, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    
    set(gca, 'XTick', 1:3, 'XTickLabel', {'View', 'Quest', 'Gray'});
    ylabel(strrep(colName, '_', ' '));
    title(titles{mi});
    xlim([0.5 3.5]);
    grid on;
    
    % 标注 n
    text(0.6, max(ylim)*0.95, sprintf('n=%d chains', size(chains,1)), 'FontSize', 9);
end

sgtitle('Three-Stage Chain: View → Questionnaire → Gray', 'FontWeight', 'bold');

figFile = fullfile(fp, sprintf('%s_three_stage_chain.png', base));
saveas(fig, figFile);
fprintf('Saved figure: %s\n', figFile);
end

function plot_scene_sequence(T, fp, base)
% Scene 序列曲线：按 scene_id 的折线图
% 一眼看出是否有场景效应/趋势/异常

conds = string(T.cond);
view_mask = conds == "view";
Tv = T(view_mask, :);

if isempty(Tv)
    warning('No view segments for scene sequence plot. Skipping.');
    return;
end

% 按 scene_id 排序
[~, ord] = sort(Tv.scene_id);
Tv = Tv(ord, :);

fig = figure('Name', 'Scene Sequence', 'Position', [100 100 1000 400]);

metrics = {'O_alpha', 'O_theta', 'F_theta', 'F_alpha'};
colors = lines(numel(metrics));

hold on;
for mi = 1:numel(metrics)
    colName = metrics{mi};
    if ~ismember(colName, Tv.Properties.VariableNames)
        continue;
    end
    
    vals = Tv.(colName);
    scene_ids = Tv.scene_id;
    
    plot(scene_ids, vals, 'o-', 'Color', colors(mi,:), 'LineWidth', 1.5, ...
        'MarkerFaceColor', colors(mi,:), 'DisplayName', strrep(colName, '_', ' '));
end

xlabel('Scene ID');
ylabel('Relative Power');
title('Power Across Scenes (View Segments)');
legend('Location', 'best');
grid on;

% 标注 block 边界
if max(Tv.block_id) > 1
    block_boundary = find(diff(Tv.block_id) > 0);
    if ~isempty(block_boundary)
        xline(Tv.scene_id(block_boundary) + 0.5, 'k--', 'LineWidth', 1.5);
        text(Tv.scene_id(block_boundary) + 0.5, max(ylim)*0.95, 'Block boundary', ...
            'FontSize', 9, 'Rotation', 90);
    end
end

figFile = fullfile(fp, sprintf('%s_scene_sequence.png', base));
saveas(fig, figFile);
fprintf('Saved figure: %s\n', figFile);
end


function cfg = load_cfg(config_path)
% 读取 JSON 配置（可选）并返回 cfg 结构体
cfg = struct();
try
    if exist(config_path, 'file')
        raw = fileread(config_path);
        cfg = jsondecode(raw);
    end
catch
    cfg = struct();
end

% 默认值
if ~isfield(cfg, 'gray_dur_min'); if ~isfield(cfg,'gray_dur_min'); cfg.gray_dur_min = 3; end end
if ~isfield(cfg, 'gray_dur_max'); if ~isfield(cfg,'gray_dur_max'); cfg.gray_dur_max = 15; end end
if ~isfield(cfg, 'quest_dur_min'); if ~isfield(cfg,'quest_dur_min'); cfg.quest_dur_min = 5; end end
if ~isfield(cfg, 'quest_dur_max'); if ~isfield(cfg,'quest_dur_max'); cfg.quest_dur_max = 120; end end
if ~isfield(cfg, 'pairing_mode'); cfg.pairing_mode = 'strict'; end
if ~isfield(cfg, 'verbose'); if ~isfield(cfg,'verbose'); cfg.verbose = true; end end
if ~isfield(cfg, 'log_file'); cfg.log_file = ''; end
end
