function fp_run = get_run_dir(fp_in, cfg)
%GET_RUN_DIR Return the active run directory under output root.
%
% New tidy layout target:
%   <output_root>/runs/run_YYYYMMDD_HHMMSS/
%
% Compatibility behavior:
% - If cfg.timestamp_output_root=true, keep using output root itself as the
%   active run dir (historical behavior).
% - Otherwise create/use:
%     <output_root>/runs/current/
%
% This helper introduces a stable abstraction layer without forcing a full
% output migration in one pass.

fp_out = pipeline.get_output_root(fp_in, cfg);

useLegacyRootAsRun = false;
try
    if isfield(cfg,'timestamp_output_root') && logical(cfg.timestamp_output_root)
        useLegacyRootAsRun = true;
    end
catch
end

if useLegacyRootAsRun
    fp_run = fp_out;
else
    fp_run = fullfile(fp_out, 'runs', 'current');
end

if ~exist(fp_run,'dir')
    mkdir(fp_run);
end
end
