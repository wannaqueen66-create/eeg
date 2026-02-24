function md_path = write_global_report_md(fp_in, cfg)
%WRITE_GLOBAL_REPORT_MD Write a markdown report for batch outputs.
%
% Creates:
%   <output_root>/summary/summary_report.md

fp_out = pipeline.get_output_root(fp_in, cfg);
fp_summary = pipeline.get_summary_dir(fp_in, cfg);
fp_rep = fp_summary;
try
    if exist('pipeline.get_report_dir','file')==2
        fp_rep = pipeline.get_report_dir(fp_summary, cfg);
    end
catch
end
md_path = fullfile(fp_rep, 'summary_report.md');

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
    d = dir(fp_out);
    d = d([d.isdir]);
    names = setdiff({d.name}, {'.','..','summary'});
    names = sort(names);
    fprintf(fid, '- n_subjects: %d\n\n', numel(names));

    for i = 1:numel(names)
        sub = names{i};
        % Try to link to subject report if exists
        rep = fullfile(fp_out, sub, sprintf('%s_report.md', sub));
        if exist(rep, 'file')
            fprintf(fid, '- %s: `%s`\n', sub, rep);
        else
            fprintf(fid, '- %s\n', sub);
        end
    end

    fprintf(fid, '\n## Global tables\n\n');
    g = fullfile(fp_summary, 'global_bandpower_summary.csv');
    if exist(g, 'file')
        fprintf(fid, '- global_bandpower_summary: `%s`\n', g);
    else
        fprintf(fid, '- global_bandpower_summary: (not generated; set `global_summary=true` to enable)\n');
    end

    fclose(fid);
catch
end

end
