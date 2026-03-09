function fp_rep = get_batch_report_dir(fp_in, cfg)
%GET_BATCH_REPORT_DIR Return batch report directory.
fp_batch = pipeline.get_batch_dir(fp_in, cfg);
fp_rep = fullfile(fp_batch, 'reports');
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
end
