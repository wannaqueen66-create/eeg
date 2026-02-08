function fp_out = resolve_output_dir(fp, output_dir)
%RESOLVE_OUTPUT_DIR Resolve relative/absolute output directory.

if startsWith(output_dir, filesep) || (~isempty(regexp(output_dir,'^[A-Za-z]:', 'once')))
    fp_out = output_dir;
else
    fp_out = fullfile(fp, output_dir);
end

end
