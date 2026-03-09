function fp_rep = get_subject_report_dir(fp_in, subject_id, cfg)
%GET_SUBJECT_REPORT_DIR Return subject-level report directory.
fp_sub = pipeline.get_subject_dir(fp_in, subject_id, cfg);
fp_rep = fullfile(fp_sub, 'report');
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
end
