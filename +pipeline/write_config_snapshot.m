function write_config_snapshot(cfg, fp)
%WRITE_CONFIG_SNAPSHOT Persist config snapshot used for this run.

try
    raw = jsonencode(cfg);
    fid = fopen(fullfile(fp, 'config_used.json'), 'w');
    if fid ~= -1
        fwrite(fid, raw);
        fclose(fid);
    end
catch
end

end
