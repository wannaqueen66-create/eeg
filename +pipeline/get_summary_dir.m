function fp_summary = get_summary_dir(fp_in, cfg)
%GET_SUMMARY_DIR Create/return the summary output directory.
%
% Summary dir lives under the output root:
%   <output_root>/summary/

fp_out = pipeline.get_output_root(fp_in, cfg);
fp_summary = fullfile(fp_out, 'summary');
if ~exist(fp_summary, 'dir')
    mkdir(fp_summary);
end

end
