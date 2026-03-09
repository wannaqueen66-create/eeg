function fp_qc = get_batch_qc_dir(fp_in, cfg)
%GET_BATCH_QC_DIR Return batch QC directory.
fp_batch = pipeline.get_batch_dir(fp_in, cfg);
fp_qc = fullfile(fp_batch, 'qc');
if ~exist(fp_qc,'dir'); mkdir(fp_qc); end
end
