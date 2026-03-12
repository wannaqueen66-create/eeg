function [fp_paper, fp_tbl, fp_fig, fp_rep] = get_paper_subdirs(fp_batch, tag)
%GET_PAPER_SUBDIRS Return staged paper-output directories.
%
% Preferred staged layout:
%   <batch>/paper/<tag>/
%     tables/
%     figures/
%   <batch>/paper/report/
%
% Legacy-compatible fallback:
%   <fp_batch>/paper_fig or <fp_batch>/paper_fig_qc

if nargin < 2 || strlength(string(tag))==0
    tag = 'raw';
end

tag = char(strtrim(string(tag)));
leaf = '';
try
    [~, leaf] = fileparts(fp_batch);
catch
end

if strcmpi(leaf,'batch')
    fp_paper = fullfile(fp_batch, 'paper', tag);
    fp_tbl = fullfile(fp_paper, 'tables');
    fp_fig = fullfile(fp_paper, 'figures');
    fp_rep = fullfile(fp_batch, 'paper', 'report');
else
    if strcmpi(tag,'qc')
        fp_paper = fullfile(fp_batch, 'paper_fig_qc');
    else
        fp_paper = fullfile(fp_batch, 'paper_fig');
    end
    fp_tbl = fp_paper;
    fp_fig = fp_paper;
    fp_rep = fp_paper;
end

if ~exist(fp_paper,'dir'); mkdir(fp_paper); end
if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
if ~exist(fp_fig,'dir'); mkdir(fp_fig); end
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
end
