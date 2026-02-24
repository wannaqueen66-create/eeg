function fp_rep = get_report_dir(fp_sum, cfg)
%GET_REPORT_DIR Return summary report directory.
% tidy:  <fp_sum>/reports/
% legacy: <fp_sum>/

layout = 'legacy';
try
    if isfield(cfg,'output_layout') && ~isempty(cfg.output_layout)
        layout = char(string(cfg.output_layout));
    end
catch
end
layout = lower(strtrim(layout));

if strcmp(layout,'tidy')
    fp_rep = fullfile(fp_sum, 'reports');
else
    fp_rep = fp_sum;
end

if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
end
