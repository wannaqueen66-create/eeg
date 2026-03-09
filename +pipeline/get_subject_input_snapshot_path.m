function fp_txt = get_subject_input_snapshot_path(fp_in, subject_id, cfg)
%GET_SUBJECT_INPUT_SNAPSHOT_PATH Return path for input_set_path.txt.
subject_id = pipeline.sanitize_filename(subject_id);
fp_sub = pipeline.get_subject_dir(fp_in, subject_id, cfg);
fp_txt = fullfile(fp_sub, 'report', 'input_set_path.txt');
end
