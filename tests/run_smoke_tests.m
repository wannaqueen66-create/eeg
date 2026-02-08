function run_smoke_tests()
%RUN_SMOKE_TESTS Lightweight smoke tests for modularized pipeline.

fprintf('== EEG pipeline smoke tests ==\n');

% Test 1: config loader defaults
cfg = pipeline.load_config('__nonexistent__.json');
assert(isfield(cfg,'gray_dur_min') && cfg.gray_dur_min==3);
assert(isfield(cfg,'pairing_mode') && strcmp(cfg.pairing_mode,'strict'));
fprintf('[OK] load_config defaults\n');

% Test 2: parse input rejects invalid extension
caught = false;
try
    pipeline.parse_input('abc.txt');
catch
    caught = true;
end
assert(caught, 'parse_input should reject non-.set files');
fprintf('[OK] parse_input invalid ext check\n');

% Test 3: output prep creates expected dirs
tmp = fullfile(tempdir, ['eeg_test_' char(java.util.UUID.randomUUID)]);
mkdir(tmp);
cleanup = onCleanup(@() safe_rmdir(tmp)); %#ok<NASGU>

cfg2 = cfg;
cfg2.output_dir = 'out';
[fp_sub, fp_csv, fp_fig, fp_qc] = pipeline.prepare_output(tmp, 'subject1', cfg2);
assert(exist(fp_sub,'dir')==7);
assert(exist(fp_csv,'dir')==7);
assert(exist(fp_fig,'dir')==7);
assert(exist(fp_qc,'dir')==7);
assert(exist(fullfile(fp_sub,'config_used.json'),'file')==2);
fprintf('[OK] prepare_output directory + config snapshot\n');

fprintf('All smoke tests passed.\n');
end

function safe_rmdir(p)
try
    if exist(p,'dir')==7
        rmdir(p,'s');
    end
catch
end
end
