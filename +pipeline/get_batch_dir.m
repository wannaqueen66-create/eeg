function fp_batch = get_batch_dir(fp_in, cfg)
%GET_BATCH_DIR Return batch-level output directory.
%
% New tidy layout target:
%   <run_dir>/batch/
%
% Compatibility behavior:
% - If cfg.timestamp_output_root=true, preserve historical summary root:
%     <output_root>/summary/
% - Otherwise use:
%     <run_dir>/batch/

useLegacyBatchLayout = false;
try
    if isfield(cfg,'timestamp_output_root') && logical(cfg.timestamp_output_root)
        useLegacyBatchLayout = true;
    end
catch
end

if useLegacyBatchLayout
    fp_batch = pipeline.get_summary_dir(fp_in, cfg);
else
    fp_run = pipeline.get_run_dir(fp_in, cfg);
    fp_batch = fullfile(fp_run, 'batch');
    if ~exist(fp_batch,'dir'); mkdir(fp_batch); end
end
end
