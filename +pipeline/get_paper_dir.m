function fp_paper = get_paper_dir(fp_in, cfg, tag)
%GET_PAPER_DIR Return paper-delivery directory.
%
% New tidy layout target:
%   <run_dir>/paper/<tag>/
%
% Compatibility behavior:
% - If cfg.timestamp_output_root=true, preserve historical summary-based layout:
%     <summary>/paper_fig or paper_fig_qc
% - Otherwise use the new paper root.

if nargin < 3; tag = ''; end
tag = char(strtrim(string(tag)));

useLegacyPaperLayout = false;
try
    if isfield(cfg,'timestamp_output_root') && logical(cfg.timestamp_output_root)
        useLegacyPaperLayout = true;
    end
catch
end

if useLegacyPaperLayout
    fp_sum = pipeline.get_summary_dir(fp_in, cfg);
    if strcmpi(tag,'qc')
        fp_paper = fullfile(fp_sum, 'paper_fig_qc');
    else
        fp_paper = fullfile(fp_sum, 'paper_fig');
    end
else
    fp_run = pipeline.get_run_dir(fp_in, cfg);
    fp_paper = fullfile(fp_run, 'paper');
    if ~isempty(tag)
        fp_paper = fullfile(fp_paper, tag);
    end
end

if ~exist(fp_paper,'dir'); mkdir(fp_paper); end
end
