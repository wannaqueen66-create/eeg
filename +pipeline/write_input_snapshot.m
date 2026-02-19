function write_input_snapshot(fp_sub, input_file)
%WRITE_INPUT_SNAPSHOT Persist the exact input .set path for traceability.
%
% Writes <fp_sub>/input_set_path.txt

try
    fid = fopen(fullfile(fp_sub, 'input_set_path.txt'), 'w');
    if fid ~= -1
        fprintf(fid, '%s\n', input_file);
        fclose(fid);
    end
catch
end

end
