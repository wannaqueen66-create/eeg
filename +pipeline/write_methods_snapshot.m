function write_methods_snapshot(fp_sum, cfg)
%WRITE_METHODS_SNAPSHOT Write a journal-friendly methods snapshot markdown.

md = fullfile(fp_sum, 'methods_snapshot.md');

try
    fid = fopen(md, 'w');
    if fid==-1
        return;
    end

    fprintf(fid, '# Methods Snapshot (Auto-generated)\n\n');
    fprintf(fid, '- Generated: %s\n', datestr(now,31));

    % Git commit (best effort)
    try
        [st, out] = system('git rev-parse HEAD');
        if st==0
            fprintf(fid, '- Git commit: `%s`\n', strtrim(out));
        end
    catch
    end

    fprintf(fid, '\n## Key Settings\n\n');
    fprintf(fid, '- pairing_mode: %s\n', string(getfield_def(cfg,'pairing_mode','strict')));
    fprintf(fid, '- strict_structure: %s\n', string(getfield_def(cfg,'strict_structure',true)));
    fprintf(fid, '- figure_visible: %s\n', string(getfield_def(cfg,'figure_visible',false)));
    fprintf(fid, '- batch_summaries: %s\n', string(getfield_def(cfg,'batch_summaries',true)));

    % Paper metrics / QC
    try
        if isfield(cfg,'paper_metrics')
            fprintf(fid, '- paper_metrics: %s\n', strjoin(string(cfg.paper_metrics), ', '));
        end
    catch
    end
    try
        if isfield(cfg,'qc_apply')
            fprintf(fid, '- qc_apply: %s\n', string(cfg.qc_apply));
            fprintf(fid, '  - qc_hf_threshold: %s\n', string(getfield_def(cfg,'qc_hf_threshold','')));
            fprintf(fid, '  - qc_view_bad_frac_threshold: %s\n', string(getfield_def(cfg,'qc_view_bad_frac_threshold','')));
            fprintf(fid, '  - qc_rms_method: %s\n', string(getfield_def(cfg,'qc_rms_method','')));
            fprintf(fid, '  - qc_rms_k: %s\n', string(getfield_def(cfg,'qc_rms_k','')));
        end
    catch
    end

    if isfield(cfg,'roi')
        fprintf(fid, '\n## ROI\n\n');
        r = cfg.roi;
        fprintf(fid, '- front: %s\n', strjoin(string(r.front), ', '));
        fprintf(fid, '- par: %s\n', strjoin(string(r.par), ', '));
        fprintf(fid, '- occ: %s\n', strjoin(string(r.occ), ', '));
    end

    if isfield(cfg,'bands')
        fprintf(fid, '\n## Bands\n\n');
        b = cfg.bands;
        keys = fieldnames(b);
        for i=1:numel(keys)
            v = b.(keys{i});
            if isnumeric(v) && numel(v)==2
                fprintf(fid, '- %s: [%g, %g] Hz\n', keys{i}, v(1), v(2));
            end
        end
    end

    fclose(fid);
catch
end

end

function v = getfield_def(s, f, d)
try
    if isfield(s,f)
        v = s.(f);
    else
        v = d;
    end
catch
    v = d;
end
end
