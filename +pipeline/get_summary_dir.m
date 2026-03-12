function fp_summary = get_summary_dir(fp_in, cfg)
%GET_SUMMARY_DIR Backward-compatible alias for the batch root directory.
%
% Historical behavior:
%   <output_root>/summary/
%
% Main-branch behavior:
% - prefer calling pipeline.get_batch_dir(...) directly in new code
% - retain this helper only as a compatibility alias for older callers

fp_summary = pipeline.get_batch_dir(fp_in, cfg);

if ~exist(fp_summary, 'dir')
    mkdir(fp_summary);
end

end
