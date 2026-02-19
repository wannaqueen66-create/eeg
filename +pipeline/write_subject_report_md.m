function md_path = write_subject_report_md(fp_sub, base, this_file, cfg, T)
%WRITE_SUBJECT_REPORT_MD Write a compact markdown report per subject.
%
% Creates:
%   <fp_sub>/<base>_report.md

md_path = fullfile(fp_sub, sprintf('%s_report.md', base));

try
    fid = fopen(md_path, 'w');
    if fid == -1
        return;
    end

    fprintf(fid, '# EEG Bandpower Report: %s\n\n', base);
    fprintf(fid, '- Generated: %s\n', datestr(now, 31));
    fprintf(fid, '- Input .set: `%s`\n', this_file);
    fprintf(fid, '- Output folder: `%s`\n\n', fp_sub);

    % Config snapshot (key items)
    fprintf(fid, '## Config (key)\n\n');
    if isfield(cfg,'pairing_mode'); fprintf(fid, '- pairing_mode: `%s`\n', string(cfg.pairing_mode)); end
    if isfield(cfg,'gray_dur_min'); fprintf(fid, '- gray_dur_min: %g\n', cfg.gray_dur_min); end
    if isfield(cfg,'gray_dur_max'); fprintf(fid, '- gray_dur_max: %g\n', cfg.gray_dur_max); end
    if isfield(cfg,'quest_dur_min'); fprintf(fid, '- quest_dur_min: %g\n', cfg.quest_dur_min); end
    if isfield(cfg,'quest_dur_max'); fprintf(fid, '- quest_dur_max: %g\n', cfg.quest_dur_max); end
    fprintf(fid, '\n');

    % Segment counts
    fprintf(fid, '## Segment counts\n\n');
    try
        conds = string(T.cond);
        uC = unique(conds, 'stable');
        for i = 1:numel(uC)
            c = uC(i);
            fprintf(fid, '- %s: %d\n', c, sum(conds==c));
        end
    catch
        fprintf(fid, '- (failed to compute segment counts)\n');
    end
    fprintf(fid, '\n');

    % Pairing summary (from exported pairs_check if exists)
    fprintf(fid, '## View–Gray pairing summary\n\n');
    pairsFile = fullfile(fp_sub, 'qc', sprintf('%s_pairs_check.csv', base));
    if exist(pairsFile, 'file')
        try
            P = readtable(pairsFile);
            fprintf(fid, '- pairs_check: `%s`\n', pairsFile);
            fprintf(fid, '- n_pairs: %d\n', height(P));
            if any(strcmp(P.Properties.VariableNames,'pair_status'))
                fprintf(fid, '- pair_status breakdown:\n');
                st = string(P.pair_status);
                uS = unique(st,'stable');
                for i=1:numel(uS)
                    fprintf(fid, '  - %s: %d\n', uS(i), sum(st==uS(i)));
                end
            end
            if any(strcmp(P.Properties.VariableNames,'gap_sec'))
                fprintf(fid, '- gap_sec mean±sd: %.2f ± %.2f\n', mean(P.gap_sec,'omitnan'), std(P.gap_sec,'omitnan'));
            end
            if any(strcmp(P.Properties.VariableNames,'delta_O_alpha'))
                fprintf(fid, '- delta_O_alpha (gray-view) mean±sd: %.4f ± %.4f\n', mean(P.delta_O_alpha,'omitnan'), std(P.delta_O_alpha,'omitnan'));
            end
        catch
            fprintf(fid, '- pairs_check exists but failed to read: `%s`\n', pairsFile);
        end
    else
        fprintf(fid, '- pairs_check not found (pairing may have been skipped).\n');
    end
    fprintf(fid, '\n');

    % Key alpha means (as already printed to console)
    fprintf(fid, '## Key alpha means (Occipital, relative)\n\n');
    try
        m_closed = mean(T.O_alpha(string(T.cond)=="eyes_closed"), 'omitnan');
        m_open   = mean(T.O_alpha(string(T.cond)=="eyes_open"),   'omitnan');
        m_view   = mean(T.O_alpha(string(T.cond)=="view"),        'omitnan');
        m_gray   = mean(T.O_alpha(string(T.cond)=="gray"),        'omitnan');
        fprintf(fid, '- eyes_closed: %.4f\n', m_closed);
        fprintf(fid, '- eyes_open: %.4f\n', m_open);
        fprintf(fid, '- view: %.4f\n', m_view);
        fprintf(fid, '- gray: %.4f\n', m_gray);
        fprintf(fid, '- gray - view: %.4f\n', m_gray - m_view);
    catch
        fprintf(fid, '- (failed to compute alpha means)\n');
    end
    fprintf(fid, '\n');

    fprintf(fid, '## Outputs\n\n');
    fprintf(fid, '- CSV: `%s`\n', fullfile(fp_sub, 'csv'));
    fprintf(fid, '- Figures: `%s`\n', fullfile(fp_sub, 'fig'));
    fprintf(fid, '- QC: `%s`\n', fullfile(fp_sub, 'qc'));

    fclose(fid);
catch
end

end
