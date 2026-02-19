function fp_out = get_output_root(fp_in, cfg)
%GET_OUTPUT_ROOT Return the root output directory for a given input folder.
%
% If cfg.output_dir is empty, default to <fp_in>/bandpower_outputs.

if isfield(cfg, 'output_dir') && ~isempty(cfg.output_dir)
    fp_out = pipeline.resolve_output_dir(fp_in, cfg.output_dir);
else
    fp_out = fullfile(fp_in, 'bandpower_outputs');
end

if ~exist(fp_out, 'dir')
    mkdir(fp_out);
end

end
