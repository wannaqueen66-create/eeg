function fp_fig = get_fig_dir(fp_sum, cfg, category, tag)
%GET_FIG_DIR Return figure output directory.
%
% Historical behavior:
%   - legacy:
%       group/recovery/topo/scene_sequence/wwr -> <fp_sum>/fig/
%       paper raw -> <fp_sum>/paper_fig/
%       paper qc  -> <fp_sum>/paper_fig_qc/
%   - tidy:
%       <fp_sum>/figures/<category>/<tag>/
%
% Stage-1 refactor behavior:
% - when fp_sum points at the new batch dir, route figures to:
%     <batch>/figures/<category>/<tag>/
% - otherwise preserve historical semantics

if nargin < 4; tag = ''; end
if nargin < 3 || strlength(string(category))==0
    category = 'group';
end

layout = 'legacy';
try
    if isfield(cfg,'output_layout') && ~isempty(cfg.output_layout)
        layout = char(string(cfg.output_layout));
    end
catch
end
layout = lower(strtrim(layout));

category = lower(strtrim(char(string(category))));
tag = lower(strtrim(char(string(tag))));

leaf = '';
try
    [~, leaf] = fileparts(fp_sum);
catch
end

if strcmp(layout,'tidy')
    if strcmpi(leaf,'batch')
        fp_fig = fullfile(fp_sum, 'figures', category);
    else
        fp_fig = fullfile(fp_sum, 'figures', category);
    end
    if ~isempty(tag)
        fp_fig = fullfile(fp_fig, tag);
    end
else
    % legacy
    if strcmp(category,'paper')
        if strcmp(tag,'qc')
            fp_fig = fullfile(fp_sum, 'paper_fig_qc');
        else
            fp_fig = fullfile(fp_sum, 'paper_fig');
        end
    else
        fp_fig = fullfile(fp_sum, 'fig');
    end
end

if ~exist(fp_fig,'dir')
    mkdir(fp_fig);
end
end
