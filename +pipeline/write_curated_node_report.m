function md_path = write_curated_node_report(nodeRoot, layerKind, branchKind)
%WRITE_CURATED_NODE_REPORT Write a small self-explanatory README for one curated node.
%
% nodeRoot example:
%   <batch>/descriptive/overall
%   <batch>/descriptive/experience
%   <batch>/inferential/overall
%   <batch>/inferential/experience
%
% The goal is lightweight readability inside each folder, so users can open
% the node directly without first reading pipeline code or repo docs.

if nargin < 3
    md_path = '';
    return;
end

md_path = '';
try
    if isempty(nodeRoot) || ~exist(nodeRoot, 'dir')
        return;
    end

    layerKind = string(layerKind);
    branchKind = string(branchKind);

    fpTables = fullfile(nodeRoot, 'tables');
    fpFigures = fullfile(nodeRoot, 'figures');
    fpReport = fullfile(nodeRoot, 'report');
    if ~exist(fpTables, 'dir'); mkdir(fpTables); end
    if ~exist(fpFigures, 'dir'); mkdir(fpFigures); end
    if ~exist(fpReport, 'dir'); mkdir(fpReport); end

    md_path = fullfile(fpReport, 'README_NODE.md');
    fid = fopen(md_path, 'w');
    if fid == -1
        md_path = '';
        return;
    end

    fprintf(fid, '# Curated Node Guide: %s / %s\n\n', char(titlecase_local(layerKind)), char(titlecase_local(branchKind)));
    fprintf(fid, '- Node root: `%s`\n', nodeRoot);
    fprintf(fid, '- Folder rule: `tables/` = CSV only; `figures/` = PNG only; `report/` = Markdown only.\n\n');

    fprintf(fid, '## What this node is for\n\n');
    switch char(layerKind) + "::" + char(branchKind)
        case "descriptive::overall"
            fprintf(fid, 'This node is the simplest full-sample descriptive entry for the redesigned main branch.\n');
            fprintf(fid, 'Use it when you want the cleanest overall means, recovery summaries, and WWR × Complexity descriptive patterns without opening task folders first.\n\n');
            fprintf(fid, 'Recommended reading order:\n');
            fprintf(fid, '1. `tables/overall_scene_metric_means_qc.csv` (or `_raw.csv` if QC version is absent)\n');
            fprintf(fid, '2. `figures/overall_metric_bar_qc.png`\n');
            fprintf(fid, '3. `tables/overall_scene_metric_means_by_WWR_Complexity_qc.csv`\n');
            fprintf(fid, '4. `figures/overall_factor_grid_qc.png`\n');
            fprintf(fid, '5. `tables/overall_recovery_means_qc.csv`\n');
            fprintf(fid, '6. `figures/overall_recovery_bar_qc.png`\n\n');
        case "descriptive::experience"
            fprintf(fid, 'This node is the default grouped descriptive entry for the redesigned main branch.\n');
            fprintf(fid, 'Use it when you want Experience High vs Low summaries, grouped recovery summaries, and Experience × WWR × Complexity descriptive patterns.\n\n');
            fprintf(fid, 'Recommended reading order:\n');
            fprintf(fid, '1. `tables/experience_scene_metric_means_qc.csv` (or `_raw.csv`)\n');
            fprintf(fid, '2. `figures/experience_metric_bar_qc.png`\n');
            fprintf(fid, '3. `tables/experience_scene_metric_means_by_WWR_Complexity_qc.csv`\n');
            fprintf(fid, '4. `figures/experience_factor_grid_qc.png`\n');
            fprintf(fid, '5. `tables/experience_recovery_means_qc.csv`\n');
            fprintf(fid, '6. `figures/experience_recovery_bar_qc.png`\n');
            fprintf(fid, '7. `tables/experience_trialindex_neural_response_qc.csv`\n');
            fprintf(fid, '8. `figures/experience_trialindex_response_qc.png`\n\n');
        case "inferential::overall"
            fprintf(fid, 'This node is the curated full-sample inferential entry for the redesigned main branch.\n');
            fprintf(fid, 'Use it when you want the main effect / interaction / trend results for the whole sample without traversing detailed task outputs first.\n\n');
            fprintf(fid, 'Recommended reading order:\n');
            fprintf(fid, '1. `tables/overall_inferential_summary_qc.csv` (or `_raw.csv`)\n');
            fprintf(fid, '2. `figures/overall_inferential_heatmap_qc.png`\n');
            fprintf(fid, '3. `tables/overall_wwr_trend_summary_qc.csv`\n');
            fprintf(fid, '4. `figures/overall_wwr_trend_heatmap_qc.png`\n');
            fprintf(fid, '5. `tables/overall_recovery_inferential_summary_qc.csv`\n');
            fprintf(fid, '6. `figures/overall_recovery_inferential_qc.png`\n\n');
        case "inferential::experience"
            fprintf(fid, 'This node is the curated experience-group inferential entry for the redesigned main branch.\n');
            fprintf(fid, 'Use it when you want Experience-group modulation, group contrasts, grouped WWR trends, and selected mirrored task outputs.\n\n');
            fprintf(fid, 'Recommended reading order:\n');
            fprintf(fid, '1. `tables/experience_inferential_summary_qc.csv` (or `_raw.csv`)\n');
            fprintf(fid, '2. `figures/experience_inferential_heatmap_qc.png`\n');
            fprintf(fid, '3. `tables/experience_wwr_trend_summary_qc.csv`\n');
            fprintf(fid, '4. `figures/experience_wwr_trend_heatmap_qc.png`\n');
            fprintf(fid, '5. `tables/experience_recovery_inferential_summary_qc.csv`\n');
            fprintf(fid, '6. `figures/experience_recovery_inferential_qc.png`\n');
            fprintf(fid, '7. `tables/trialindex_lmm_summary_qc.csv` (mirrored from task3 when available)\n');
            fprintf(fid, '8. `figures/trialindex_lmm_overview_qc_experience.png` (mirrored when available)\n');
            fprintf(fid, '9. `tables/experience_inferential_file_index.csv` for detailed source tracing\n\n');
        otherwise
            fprintf(fid, 'This curated node provides a simplified reading surface for the redesigned main branch.\n\n');
    end

    fprintf(fid, '## Current contents snapshot\n\n');
    list_dir_md(fid, fpTables, 'CSV tables');
    list_dir_md(fid, fpFigures, 'PNG figures');
    list_dir_md(fid, fpReport, 'Markdown reports');

    fprintf(fid, '\n## Notes\n\n');
    fprintf(fid, '- Prefer QC files when both `_raw` and `_qc` exist and the use case is formal reporting.\n');
    fprintf(fid, '- `SportFreq` is intentionally not promoted to a first-class curated branch at this stage; look in detailed `analysis/` outputs if needed.\n');
    fprintf(fid, '- Detailed task-oriented outputs remain under `batch/analysis/`. This curated node is only the cleaner front surface.\n');

    fclose(fid);
catch
    try
        fclose(fid);
    catch
    end
    md_path = '';
end

end

function s = titlecase_local(x)
s = string(x);
if strlength(s) == 0
    return;
end
c = char(lower(s));
c(1) = upper(c(1));
s = string(c);
end

function list_dir_md(fid, folderPath, label)
fprintf(fid, '### %s\n\n', label);
if ~exist(folderPath, 'dir')
    fprintf(fid, '- (folder missing)\n\n');
    return;
end
D = dir(folderPath);
D = D(~[D.isdir]);
if isempty(D)
    fprintf(fid, '- (empty)\n\n');
    return;
end
[~, ord] = sort(lower({D.name}));
D = D(ord);
for i = 1:numel(D)
    fprintf(fid, '- `%s`\n', D(i).name);
end
fprintf(fid, '\n');
end
