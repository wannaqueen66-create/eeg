function fp_tbl = get_table_dir(fp_sum, cfg, category)
%GET_TABLE_DIR Return table directory for a batch-level category.
%
% Historical behavior:
%   tidy:  <fp_sum>/tables/<category>/
%   legacy: <fp_sum>/
%
% Stage-1 refactor behavior:
% - when fp_sum points at the new batch dir, route common categories to
%   batch/merged or batch/qc directly
% - otherwise preserve historical semantics

if nargin < 3 || strlength(string(category))==0
    category = 'merged_raw';
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

leaf = '';
try
    [~, leaf] = fileparts(fp_sum);
catch
end

if strcmp(layout,'tidy') && strcmpi(leaf,'batch')
    switch category
        case 'merged_raw'
            fp_tbl = fullfile(fp_sum, 'merged');
        case 'merged_qc'
            fp_tbl = fullfile(fp_sum, 'qc');
        case 'sensitivity'
            fp_tbl = fullfile(fp_sum, 'reports');
        otherwise
            fp_tbl = fullfile(fp_sum, 'tables', category);
    end
elseif strcmp(layout,'tidy')
    fp_tbl = fullfile(fp_sum, 'tables', category);
else
    fp_tbl = fp_sum;
end

if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
end
