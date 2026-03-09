function [fp_sub, fp_csv, fp_fig, fp_qc] = prepare_output(fp, subject_id, cfg)
%PREPARE_OUTPUT Create output directories and write config snapshot.
%
% Stage-1 refactor behavior:
% - preserve historical behavior when timestamp_output_root=true
% - otherwise route subject outputs through the new subject path helpers
%   so later directory migrations can happen centrally.

fp_sub = pipeline.get_subject_dir(fp, subject_id, cfg);
fp_csv = pipeline.get_subject_table_dir(fp, subject_id, cfg);
fp_fig = pipeline.get_subject_figure_dir(fp, subject_id, cfg);
fp_qc  = pipeline.get_subject_qc_dir(fp, subject_id, cfg);

pipeline.write_config_snapshot(cfg, fp_sub);

end
