function md_path = write_global_report_md(fp_in, cfg)
%WRITE_GLOBAL_REPORT_MD Write a markdown report for batch outputs.
%
% Creates:
%   staged layout: <batch>/reports/batch_report.md
%   legacy-compatible alias: <...>/summary_report.md

fp_out = pipeline.get_output_root(fp_in, cfg);
fp_batch = pipeline.get_summary_dir(fp_in, cfg);
fp_rep = fp_batch;
try
    if exist('pipeline.get_batch_report_dir','file')==2
        fp_rep = pipeline.get_batch_report_dir(fp_in, cfg);
    elseif exist('pipeline.get_report_dir','file')==2
        fp_rep = pipeline.get_report_dir(fp_batch, cfg);
    end
catch
end
md_path = fullfile(fp_rep, 'batch_report.md');
fp_md_legacy = fullfile(fp_rep, 'summary_report.md');

try
    fid = fopen(md_path, 'w');
    if fid == -1
        return;
    end

    fprintf(fid, '# EEG Bandpower Batch Summary\n\n');
    fprintf(fid, '- Generated: %s\n', datestr(now, 31));
    fprintf(fid, '- Input folder: `%s`\n', fp_in);
    fprintf(fid, '- Output root: `%s`\n\n', fp_out);

    fprintf(fid, '## Subjects\n\n');
    names = cellstr(pipeline.list_subject_ids(fp_in, cfg));
    fprintf(fid, '- n_subjects: %d\n\n', numel(names));

    for i = 1:numel(names)
        sub = names{i};
        rep = fullfile(pipeline.get_subject_report_dir(fp_in, sub, cfg), sprintf('%s_report.md', sub));
        if ~exist(rep, 'file')
            rep = fullfile(pipeline.get_subject_dir(fp_in, sub, cfg), sprintf('%s_report.md', sub));
        end
        if exist(rep, 'file')
            fprintf(fid, '- %s: `%s`\n', sub, rep);
        else
            fprintf(fid, '- %s\n', sub);
        end
    end

    fprintf(fid, '\n## Global tables\n\n');
    g = fullfile(pipeline.get_batch_merged_dir(fp_in, cfg), 'global_bandpower_summary.csv');
    if ~exist(g, 'file')
        g = fullfile(fp_batch, 'global_bandpower_summary.csv');
    end
    if exist(g, 'file')
        fprintf(fid, '- global_bandpower_summary: `%s`\n', g);
    else
        fprintf(fid, '- global_bandpower_summary: (not generated; set `global_summary=true` to enable)\n');
    end

    fclose(fid);

    % Backward-compatible alias for older tooling/docs that still expect summary_report.md
    try
        if exist(md_path, 'file')
            copyfile(md_path, fp_md_legacy);
        end
    catch
    end
catch
end

end
