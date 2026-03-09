function fp_tbl = get_batch_merged_dir(fp_in, cfg)
%GET_BATCH_MERGED_DIR Return batch merged-table directory.
fp_batch = pipeline.get_batch_dir(fp_in, cfg);
fp_tbl = fullfile(fp_batch, 'merged');
if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
end
