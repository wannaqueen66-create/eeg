function fp_log = get_log_dir(fp_batch, cfg)
%GET_LOG_DIR Return summary logs directory.
% tidy:  <fp_batch>/logs/
% legacy: <fp_batch>/

layout = 'legacy';
try
    if isfield(cfg,'output_layout') && ~isempty(cfg.output_layout)
        layout = char(string(cfg.output_layout));
    end
catch
end
layout = lower(strtrim(layout));

if strcmp(layout,'tidy')
    fp_log = fullfile(fp_batch, 'logs');
else
    fp_log = fp_batch;
end

if ~exist(fp_log,'dir'); mkdir(fp_log); end
end
