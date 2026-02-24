function fp_aud = get_audit_dir(fp_sum, cfg)
%GET_AUDIT_DIR Return summary audit directory.
% tidy:  <fp_sum>/audit/
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
    fp_aud = fullfile(fp_sum, 'audit');
else
    fp_aud = fp_sum;
end

if ~exist(fp_aud,'dir'); mkdir(fp_aud); end
end
