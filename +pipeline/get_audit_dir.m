function fp_aud = get_audit_dir(fp_batch, cfg)
%GET_AUDIT_DIR Return summary audit directory.
% tidy:  <fp_batch>/audit/
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
    fp_aud = fullfile(fp_batch, 'audit');
else
    fp_aud = fp_batch;
end

if ~exist(fp_aud,'dir'); mkdir(fp_aud); end
end
