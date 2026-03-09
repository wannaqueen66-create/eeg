function [fp_task, fp_tbl, fp_fig, fp_rep] = get_analysis_task_subdirs(fp_sum, taskName, tag)
%GET_ANALYSIS_TASK_SUBDIRS Return staged analysis task subdirectories.
%
% Preferred staged layout:
%   <batch>/analysis/<taskName>/<tag>/
%     tables/
%     figures/
%     reports/
%
% Legacy-compatible fallback:
%   <fp_sum>/analysis-2/<taskName>/
%     tables/<tag>/
%     figures/<tag>/
%     reports/<tag>/

if nargin < 3 || strlength(string(tag))==0
    tag = 'raw';
end
if nargin < 2 || strlength(string(taskName))==0
    taskName = 'task_unknown';
end

tag = char(strtrim(string(tag)));
taskName = char(strtrim(string(taskName)));

leaf = '';
try
    [~, leaf] = fileparts(fp_sum);
catch
end

if strcmpi(leaf,'batch')
    fp_task = fullfile(fp_sum, 'analysis', taskName, tag);
    fp_tbl = fullfile(fp_task, 'tables');
    fp_fig = fullfile(fp_task, 'figures');
    fp_rep = fullfile(fp_task, 'reports');
else
    fp_task = fullfile(fp_sum, 'analysis-2', taskName);
    fp_tbl = fullfile(fp_task, 'tables', tag);
    fp_fig = fullfile(fp_task, 'figures', tag);
    fp_rep = fullfile(fp_task, 'reports', tag);
end

if ~exist(fp_task,'dir'); mkdir(fp_task); end
if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
if ~exist(fp_fig,'dir'); mkdir(fp_fig); end
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
end
