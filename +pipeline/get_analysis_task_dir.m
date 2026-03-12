function fp_task = get_analysis_task_dir(fp_in, cfg, taskName, tag)
%GET_ANALYSIS_TASK_DIR Return analysis task directory.
%
% New tidy layout target:
%   <run_dir>/analysis/<taskName>/<tag>/
%
% Compatibility behavior:
% - If cfg.timestamp_output_root=true, preserve historical summary-based layout:
%     <summary>/analysis-2/<taskName>/
%     and append <tag> only when non-empty under the caller's own subfolders.
% - Otherwise use the new analysis root.

if nargin < 4; tag = ''; end
if nargin < 3 || strlength(string(taskName))==0
    taskName = 'task_unknown';
end

taskName = char(strtrim(string(taskName)));
tag = char(strtrim(string(tag)));

useLegacyAnalysisLayout = false;
try
    if isfield(cfg,'timestamp_output_root') && logical(cfg.timestamp_output_root)
        useLegacyAnalysisLayout = true;
    end
catch
end

if useLegacyAnalysisLayout
    fp_batch = pipeline.get_batch_dir(fp_in, cfg);
    fp_task = fullfile(fp_batch, 'analysis-2', taskName);
else
    fp_run = pipeline.get_run_dir(fp_in, cfg);
    fp_task = fullfile(fp_run, 'analysis', taskName);
    if ~isempty(tag)
        fp_task = fullfile(fp_task, tag);
    end
end

if ~exist(fp_task,'dir'); mkdir(fp_task); end
end
