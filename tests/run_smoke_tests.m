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
assert(exist(fullfile(fp_sub,'report'),'dir')==7 || exist(fullfile(fp_sub,'config_used.json'),'file')==2);
assert(exist(fullfile(fp_sub,'report','config_used.json'),'file')==2 || exist(fullfile(fp_sub,'config_used.json'),'file')==2);
fprintf('[OK] prepare_output directory + config snapshot\n');

% Test 3b: staged helper directories exist
fp_run = pipeline.get_run_dir(tmp, cfg2);
fp_sub2 = pipeline.get_subject_dir(tmp, 'subject1', cfg2);
fp_tbl2 = pipeline.get_subject_table_dir(tmp, 'subject1', cfg2);
fp_fig2 = pipeline.get_subject_figure_dir(tmp, 'subject1', cfg2);
fp_qc2  = pipeline.get_subject_qc_dir(tmp, 'subject1', cfg2);
fp_rep2 = pipeline.get_subject_report_dir(tmp, 'subject1', cfg2);
fp_batch = pipeline.get_batch_dir(tmp, cfg2);
fp_bm = pipeline.get_batch_merged_dir(tmp, cfg2);
fp_bq = pipeline.get_batch_qc_dir(tmp, cfg2);
fp_br = pipeline.get_batch_report_dir(tmp, cfg2);
fp_ba = pipeline.get_batch_audit_dir(tmp, cfg2);
assert(exist(fp_run,'dir')==7);
assert(exist(fp_sub2,'dir')==7);
assert(exist(fp_tbl2,'dir')==7);
assert(exist(fp_fig2,'dir')==7);
assert(exist(fp_qc2,'dir')==7);
assert(exist(fp_rep2,'dir')==7);
assert(exist(fp_batch,'dir')==7);
assert(exist(fp_bm,'dir')==7);
assert(exist(fp_bq,'dir')==7);
assert(exist(fp_br,'dir')==7);
assert(exist(fp_ba,'dir')==7);
fprintf('[OK] staged output helper directories\n');


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

% Test 6: curated node report helper writes README_NODE.md into report/
nodeRoot = fullfile(fp_batch, 'descriptive', 'overall');
mkdir(nodeRoot);
mdNode = pipeline.write_curated_node_report(nodeRoot, 'descriptive', 'overall');
assert(~isempty(mdNode) && exist(mdNode,'file')==2, 'README_NODE.md should be generated');
nodeText = fileread(mdNode);
assert(contains(nodeText, 'Curated Node Guide: Descriptive / Overall'), 'Node guide title missing');
assert(contains(nodeText, 'Folder rule'), 'Node guide folder rule missing');
assert(contains(nodeText, 'Recommended reading order'), 'Node guide reading order missing');
fprintf('[OK] curated node report writer\n');

% Test 7: curated builder remains structurally available, but is optional in current main
curatedBuilder = fileread(fullfile(fileparts(mfilename('fullpath')), '..', '+pipeline', 'build_curated_main_outputs.m'));
assert(contains(curatedBuilder, "write_curated_node_report(fp_do, 'descriptive', 'overall')"), 'Missing descriptive/overall node report call');
assert(contains(curatedBuilder, "write_curated_node_report(fp_de, 'descriptive', 'experience')"), 'Missing descriptive/experience node report call');
assert(contains(curatedBuilder, "write_curated_node_report(fp_io, 'inferential', 'overall')"), 'Missing inferential/overall node report call');
assert(contains(curatedBuilder, "write_curated_node_report(fp_ie, 'inferential', 'experience')"), 'Missing inferential/experience node report call');
assert(contains(mainScript, 'enable_curated_main_outputs') || contains(fileread(fullfile(fileparts(mfilename('fullpath')), '..', '+pipeline', 'load_config.m')), 'enable_curated_main_outputs'), 'Expected curated output gating in current main');
fprintf('[OK] curated main builder available + gated by config\n');

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
