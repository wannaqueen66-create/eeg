function write_curated_readme_index(fp_sum)
%WRITE_CURATED_README_INDEX Write a lightweight curated entry note under batch/reports.
%
% This does not replace detailed task reports. It gives users a stable first
% page that points them to the redesigned main output surface.

if nargin < 1 || isempty(fp_sum) || ~exist(fp_sum,'dir')
    return;
end

fp_rep = fullfile(fp_sum, 'reports');
if ~exist(fp_rep,'dir'); mkdir(fp_rep); end
fp_md = fullfile(fp_rep, 'README_MAIN_ENTRY.md');

try
    fid = fopen(fp_md, 'w');
    if fid < 0
        return;
    end
    fprintf(fid, '# Main Output Entry\n\n');
    fprintf(fid, 'Recommended reading order for the redesigned main branch:\n\n');
    fprintf(fid, '1. `../descriptive/overall/`\n');
    fprintf(fid, '2. `../descriptive/experience/`\n');
    fprintf(fid, '3. `../inferential/overall/`\n');
    fprintf(fid, '4. `../inferential/experience/`\n\n');
    fprintf(fid, 'Detailed task outputs still exist under `../analysis/` for deep inspection, but the default curated surface now prioritizes the four folders above.\n');
    fclose(fid);
catch
end
end
