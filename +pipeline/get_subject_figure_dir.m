function fp_fig = get_subject_figure_dir(fp_in, subject_id, cfg)
%GET_SUBJECT_FIGURE_DIR Return subject-level figures directory.
fp_sub = pipeline.get_subject_dir(fp_in, subject_id, cfg);
fp_fig = fullfile(fp_sub, 'figures');
if ~exist(fp_fig,'dir'); mkdir(fp_fig); end
end
