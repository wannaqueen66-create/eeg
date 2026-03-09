function fp_aud = get_batch_audit_dir(fp_in, cfg)
%GET_BATCH_AUDIT_DIR Return batch audit directory.
fp_batch = pipeline.get_batch_dir(fp_in, cfg);
fp_aud = fullfile(fp_batch, 'audit');
if ~exist(fp_aud,'dir'); mkdir(fp_aud); end
end
