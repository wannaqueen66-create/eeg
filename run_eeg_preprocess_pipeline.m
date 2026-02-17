function run_eeg_preprocess_pipeline(input_path, config_path)
%RUN_EEG_PREPROCESS_PIPELINE Preprocess EEG .set files for downstream bandpower analysis.
%
% Usage:
%   run_eeg_preprocess_pipeline();
%   run_eeg_preprocess_pipeline('path/to/file_or_folder');
%   run_eeg_preprocess_pipeline('path/to/file_or_folder', 'preprocess_config.json');
%
% This script performs (configurable):
%   1) Bandpass filter
%   2) Notch filter (line noise)
%   3) Re-reference
%   4) ICA (optional)
%   5) ICLabel-based artifact IC rejection (optional, if plugin available)
%   6) Save *_preproc.set

if nargin < 1, input_path = ''; end
if nargin < 2 || isempty(config_path)
    config_path = fullfile(fileparts(mfilename('fullpath')), 'preprocess_config.json');
end

cfg = load_preproc_config(config_path);
files = pipeline.parse_input(input_path);

fprintf('== EEG preprocess pipeline ==\n');
fprintf('Files to process: %d\n', numel(files));

for i = 1:numel(files)
    inFile = files{i};
    [fp, base] = fileparts(inFile);
    fprintf('\n[%d/%d] %s\n', i, numel(files), inFile);

    stage = 'load';
    try
        EEG = pop_loadset('filename', [base '.set'], 'filepath', fp);
        EEG = eeg_checkset(EEG);

        stage = 'prep_output';
        outRoot = resolve_preproc_output(fp, cfg.output_dir);
        if exist(outRoot, 'dir') ~= 7
            mkdir(outRoot);
        end

        stage = 'bandpass';
        if cfg.do_bandpass
            EEG = apply_bandpass(EEG, cfg);
        end

        stage = 'notch';
        if cfg.do_notch
            EEG = apply_notch(EEG, cfg);
        end

        stage = 'reref';
        if cfg.do_reref
            EEG = apply_reref(EEG, cfg);
        end

        stage = 'ica';
        if cfg.do_ica
            EEG = apply_ica(EEG, cfg);
        end

        stage = 'iclabel';
        if cfg.do_iclabel_reject
            EEG = apply_iclabel_reject(EEG, cfg);
        end

        stage = 'save';
        outSet = sprintf('%s_preproc.set', base);
        EEG = pop_saveset(EEG, 'filename', outSet, 'filepath', outRoot);

        stage = 'report';
        write_preproc_report(EEG, fullfile(outRoot, sprintf('%s_preproc_report.csv', base)), cfg, inFile, fullfile(outRoot, outSet));

        fprintf('[OK] Saved: %s\n', fullfile(outRoot, outSet));
    catch ME
        fprintf('[FAIL] %s (stage=%s): %s\n', inFile, stage, ME.message);
    end
end

fprintf('\nPreprocessing finished.\n');
end

function cfg = load_preproc_config(config_path)
cfg = struct();
try
    if exist(config_path, 'file') == 2
        cfg = jsondecode(fileread(config_path));
    end
catch
    cfg = struct();
end

% defaults
if ~isfield(cfg, 'output_dir'); cfg.output_dir = 'preprocessed'; end
if ~isfield(cfg, 'do_bandpass'); cfg.do_bandpass = true; end
if ~isfield(cfg, 'bandpass_low'); cfg.bandpass_low = 1; end
if ~isfield(cfg, 'bandpass_high'); cfg.bandpass_high = 45; end
if ~isfield(cfg, 'do_notch'); cfg.do_notch = true; end
if ~isfield(cfg, 'line_freq'); cfg.line_freq = 50; end
if ~isfield(cfg, 'notch_bw'); cfg.notch_bw = 2; end
if ~isfield(cfg, 'do_reref'); cfg.do_reref = true; end
if ~isfield(cfg, 'reref_mode'); cfg.reref_mode = 'average'; end % average|channels
if ~isfield(cfg, 'reref_channels'); cfg.reref_channels = []; end
if ~isfield(cfg, 'do_ica'); cfg.do_ica = true; end
if ~isfield(cfg, 'ica_extended'); cfg.ica_extended = 1; end
if ~isfield(cfg, 'ica_pca'); cfg.ica_pca = 0; end % 0 means no PCA reduction
if ~isfield(cfg, 'do_iclabel_reject'); cfg.do_iclabel_reject = false; end
if ~isfield(cfg, 'iclabel_eye_thr'); cfg.iclabel_eye_thr = 0.9; end
if ~isfield(cfg, 'iclabel_muscle_thr'); cfg.iclabel_muscle_thr = 0.9; end
if ~isfield(cfg, 'iclabel_heart_thr'); cfg.iclabel_heart_thr = 0.9; end
end

function out = resolve_preproc_output(fp, output_dir)
if isempty(output_dir)
    out = fp;
    return;
end
if startsWith(output_dir, filesep) || (~isempty(regexp(output_dir, '^[A-Za-z]:', 'once')))
    out = output_dir;
else
    out = fullfile(fp, output_dir);
end
end

function EEG = apply_bandpass(EEG, cfg)
if exist('pop_eegfiltnew', 'file') ~= 2
    error('pop_eegfiltnew not found. Please ensure EEGLAB is on path.');
end
EEG = pop_eegfiltnew(EEG, cfg.bandpass_low, cfg.bandpass_high);
EEG = eeg_checkset(EEG);
fprintf('  - Bandpass: %.2f-%.2f Hz\n', cfg.bandpass_low, cfg.bandpass_high);
end

function EEG = apply_notch(EEG, cfg)
if exist('pop_eegfiltnew', 'file') ~= 2
    error('pop_eegfiltnew not found. Please ensure EEGLAB is on path.');
end
lo = cfg.line_freq - cfg.notch_bw/2;
hi = cfg.line_freq + cfg.notch_bw/2;
EEG = pop_eegfiltnew(EEG, lo, hi, [], 1); % revfilt=1 notch/band-stop
EEG = eeg_checkset(EEG);
fprintf('  - Notch: %.2f-%.2f Hz\n', lo, hi);
end

function EEG = apply_reref(EEG, cfg)
if strcmpi(cfg.reref_mode, 'average')
    EEG = pop_reref(EEG, []);
    fprintf('  - Re-reference: average\n');
elseif strcmpi(cfg.reref_mode, 'channels')
    if isempty(cfg.reref_channels)
        error('reref_mode=channels but reref_channels is empty');
    end
    EEG = pop_reref(EEG, cfg.reref_channels);
    fprintf('  - Re-reference: channels [%s]\n', num2str(cfg.reref_channels));
else
    error('Unsupported reref_mode: %s', cfg.reref_mode);
end
EEG = eeg_checkset(EEG);
end

function EEG = apply_ica(EEG, cfg)
if exist('pop_runica', 'file') ~= 2
    error('pop_runica not found. Please ensure EEGLAB is on path.');
end
if cfg.ica_pca > 0
    EEG = pop_runica(EEG, 'extended', cfg.ica_extended, 'pca', cfg.ica_pca, 'interrupt', 'off');
    fprintf('  - ICA: runica (extended=%d, pca=%d)\n', cfg.ica_extended, cfg.ica_pca);
else
    EEG = pop_runica(EEG, 'extended', cfg.ica_extended, 'interrupt', 'off');
    fprintf('  - ICA: runica (extended=%d)\n', cfg.ica_extended);
end
EEG = eeg_checkset(EEG);
end

function EEG = apply_iclabel_reject(EEG, cfg)
if ~isfield(EEG, 'icaweights') || isempty(EEG.icaweights)
    warning('ICLabel reject skipped: ICA not available.');
    return;
end

if exist('pop_iclabel', 'file') ~= 2
    warning('ICLabel plugin not found. Skip auto IC rejection.');
    return;
end

EEG = pop_iclabel(EEG, 'default');
cls = EEG.etc.ic_classification.ICLabel.classifications;
% ICLabel order: brain, muscle, eye, heart, line noise, channel noise, other
rej = cls(:,2) >= cfg.iclabel_muscle_thr | ...
      cls(:,3) >= cfg.iclabel_eye_thr    | ...
      cls(:,4) >= cfg.iclabel_heart_thr;
rejIdx = find(rej);

if ~isempty(rejIdx)
    EEG = pop_subcomp(EEG, rejIdx, 0);
    fprintf('  - ICLabel reject: removed %d ICs\n', numel(rejIdx));
else
    fprintf('  - ICLabel reject: removed 0 ICs\n');
end
EEG = eeg_checkset(EEG);
end

function write_preproc_report(EEG, outCsv, cfg, inFile, outFile)
subject = string('');
if isfield(EEG, 'subject') && ~isempty(EEG.subject)
    subject = string(EEG.subject);
end

nIC = 0;
if isfield(EEG, 'icaweights') && ~isempty(EEG.icaweights)
    nIC = size(EEG.icaweights, 1);
end

T = table( ...
    string(inFile), string(outFile), subject, EEG.nbchan, EEG.pnts, EEG.srate, EEG.trials, ...
    cfg.do_bandpass, cfg.bandpass_low, cfg.bandpass_high, ...
    cfg.do_notch, cfg.line_freq, cfg.notch_bw, ...
    cfg.do_reref, string(cfg.reref_mode), ...
    cfg.do_ica, nIC, ...
    cfg.do_iclabel_reject, ...
    'VariableNames', { ...
    'input_file','output_file','subject','nbchan','pnts','srate','trials', ...
    'do_bandpass','bandpass_low','bandpass_high', ...
    'do_notch','line_freq','notch_bw', ...
    'do_reref','reref_mode', ...
    'do_ica','n_ica_components', ...
    'do_iclabel_reject'});

writetable(T, outCsv);
end
