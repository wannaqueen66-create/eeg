function fp_summary = get_summary_dir(fp_in, cfg)
%GET_SUMMARY_DIR Create/return the batch-summary root directory.
%
% Historical behavior:
%   <output_root>/summary/
%
% Stage-1 refactor behavior:
% - preserve historical summary root when timestamp_output_root=true
% - otherwise route summary semantics to the new batch directory

useLegacySummary = false;
try
    if isfield(cfg,'timestamp_output_root') && logical(cfg.timestamp_output_root)
        useLegacySummary = true;
    end
catch
end

if useLegacySummary
    fp_out = pipeline.get_output_root(fp_in, cfg);
    fp_summary = fullfile(fp_out, 'summary');
else
    fp_summary = pipeline.get_batch_dir(fp_in, cfg);
end

if ~exist(fp_summary, 'dir')
    mkdir(fp_summary);
end

end
