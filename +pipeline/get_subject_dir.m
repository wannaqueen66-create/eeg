function fp_sub = get_subject_dir(fp_in, subject_id, cfg)
%GET_SUBJECT_DIR Return subject root directory for single-subject outputs.
%
% New tidy layout target:
%   <run_dir>/subjects/<subject_id>/
%
% Compatibility behavior:
% - If cfg.timestamp_output_root=true, preserve historical layout:
%     <output_root>/<subject_id>/
% - Otherwise use the new subject-oriented structure under runs/current.

subject_id = pipeline.sanitize_filename(subject_id);
fp_out = pipeline.get_output_root(fp_in, cfg);

useLegacySubjectLayout = false;
try
    if isfield(cfg,'timestamp_output_root') && logical(cfg.timestamp_output_root)
        useLegacySubjectLayout = true;
    end
catch
end

if useLegacySubjectLayout
    fp_sub = fullfile(fp_out, subject_id);
else
    fp_run = pipeline.get_run_dir(fp_in, cfg);
    fp_sub = fullfile(fp_run, 'subjects', subject_id);
end

if ~exist(fp_sub,'dir')
    mkdir(fp_sub);
end
end
