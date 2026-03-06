function run_smoke_tests()
%RUN_SMOKE_TESTS Lightweight smoke tests for modularized pipeline.

fprintf('== EEG pipeline smoke tests ==\n');

% Test 1: config loader defaults
cfg = pipeline.load_config('__nonexistent__.json');
assert(isfield(cfg,'gray_dur_min') && cfg.gray_dur_min==3);
assert(isfield(cfg,'pairing_mode') && strcmp(cfg.pairing_mode,'strict'));
assert(isfield(cfg,'eye_merge_enabled') && cfg.eye_merge_enabled==false);
assert(isfield(cfg,'eye_summary_path') && strcmp(cfg.eye_summary_path,''));
assert(isfield(cfg,'eye_merge_keys'));
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


% Test 4: scene_level export schema fields are present in pipeline script (static smoke check)
mainScript = fileread(fullfile(fileparts(mfilename('fullpath')), '..', 'run_eeg_bandpower_pipeline.m'));
requiredFields = {
    "'block_id'", "'cycle_in_block'", "'pair_id'", ...
    "'F_low_beta'", "'F_high_beta'", "'F_low_gamma'", ...
    "'P_low_beta'", "'P_high_beta'", "'P_low_gamma'", ...
    "'O_low_beta'", "'O_high_beta'", "'O_low_gamma'", ...
    "'F_TAR'", "'F_TBR'", "'F_BA'", ...
    "'P_TAR'", "'P_TBR'", "'P_BA'", ...
    "'O_TAR'", "'O_TBR'", "'O_BA'"};
for i = 1:numel(requiredFields)
    assert(contains(mainScript, requiredFields{i}), ['Missing scene_level field marker: ' requiredFields{i}]);
end
assert(contains(mainScript, 'bands.low_gamma'), 'low_gamma band support is expected in main pipeline');
fprintf('[OK] scene_level schema + low_gamma static checks\n');

% Test 5: low_gamma should be exported from out_rel40 (not out_rel30)
assert(contains(mainScript, 'T.F_low_gamma = out_rel40(:,6);'), 'F_low_gamma should use out_rel40 denominator');
assert(contains(mainScript, 'T.P_low_gamma = out_rel40(:,12);'), 'P_low_gamma should use out_rel40 denominator');
assert(contains(mainScript, 'T.O_low_gamma = out_rel40(:,18);'), 'O_low_gamma should use out_rel40 denominator');
fprintf('[OK] low_gamma denominator mapping checks\n');

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
