function fp_tbl = get_table_dir(fp_sum, cfg, category)
%GET_TABLE_DIR Return summary table directory.
% tidy:  <fp_sum>/tables/<category>/
% legacy: <fp_sum>/

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

if strcmp(layout,'tidy')
    fp_tbl = fullfile(fp_sum, 'tables', category);
else
    fp_tbl = fp_sum;
end

if ~exist(fp_tbl,'dir'); mkdir(fp_tbl); end
end
