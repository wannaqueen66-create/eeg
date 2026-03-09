function fp_json = get_subject_config_snapshot_path(fp_in, subject_id, cfg)
%GET_SUBJECT_CONFIG_SNAPSHOT_PATH Return path for config_used.json.
subject_id = pipeline.sanitize_filename(subject_id);
fp_sub = pipeline.get_subject_dir(fp_in, subject_id, cfg);
fp_json = fullfile(fp_sub, 'report', 'config_used.json');
end
