function fp_qc = get_subject_qc_dir(fp_in, subject_id, cfg)
%GET_SUBJECT_QC_DIR Return subject-level QC directory.
fp_sub = pipeline.get_subject_dir(fp_in, subject_id, cfg);
fp_qc = fullfile(fp_sub, 'qc');
if ~exist(fp_qc,'dir'); mkdir(fp_qc); end
end
