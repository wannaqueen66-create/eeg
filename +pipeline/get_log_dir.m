function fp_log = get_log_dir(fp_sum, cfg)
%GET_LOG_DIR Return summary logs directory.
% tidy:  <fp_sum>/logs/
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
    fp_log = fullfile(fp_sum, 'logs');
else
    fp_log = fp_sum;
end

if ~exist(fp_log,'dir'); mkdir(fp_log); end
end
