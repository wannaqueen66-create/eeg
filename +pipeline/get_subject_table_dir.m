function fp_tbl = get_subject_table_dir(fp_in, subject_id, cfg)
%GET_SUBJECT_TABLE_DIR Return subject-level tables directory.
fp_sub = pipeline.get_subject_dir(fp_in, subject_id, cfg);
fp_tbl = fullfile(fp_sub, 'tables');
if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
end
