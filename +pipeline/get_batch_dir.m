function fp_batch = get_batch_dir(fp_in, cfg)
%GET_BATCH_DIR Return batch-level output directory.
%
% New tidy layout target:
%   <run_dir>/batch/
%
% Compatibility behavior:
% - If cfg.timestamp_output_root=true, preserve historical summary-style root
%   as a backward-compatible fallback:
%     <output_root>/summary/
% - Otherwise use the canonical main-branch path:
%     <run_dir>/batch/

useLegacyBatchLayout = false;
try
    if isfield(cfg,'timestamp_output_root') && logical(cfg.timestamp_output_root)
        useLegacyBatchLayout = true;
    end
catch
end

if useLegacyBatchLayout
    fp_out = pipeline.get_output_root(fp_in, cfg);
    fp_batch = fullfile(fp_out, 'summary');
else
    fp_run = pipeline.get_run_dir(fp_in, cfg);
    fp_batch = fullfile(fp_run, 'batch');
end

if ~exist(fp_batch,'dir'); mkdir(fp_batch); end
end
