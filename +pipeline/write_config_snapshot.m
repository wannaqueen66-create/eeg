function write_config_snapshot(cfg, fp)
%WRITE_CONFIG_SNAPSHOT Persist config snapshot used for this run.
%
% Staged layout prefers:
%   <fp>/report/config_used.json
% Legacy fallback:
%   <fp>/config_used.json

try
    raw = jsonencode(cfg);
    fp_json = fullfile(fp, 'config_used.json');
    fp_rep = fullfile(fp, 'report');
    if exist(fp_rep,'dir')
        fp_json = fullfile(fp_rep, 'config_used.json');
    end
    fid = fopen(fp_json, 'w');
    if fid ~= -1
        fwrite(fid, raw);
        fclose(fid);
    end
catch
end

end
