function write_input_snapshot(fp_sub, input_file)
%WRITE_INPUT_SNAPSHOT Persist the exact input .set path for traceability.
%
% Staged layout prefers:
%   <fp_sub>/report/input_set_path.txt
% Legacy fallback:
%   <fp_sub>/input_set_path.txt

try
    fp_txt = fullfile(fp_sub, 'input_set_path.txt');
    fp_rep = fullfile(fp_sub, 'report');
    if exist(fp_rep,'dir')
        fp_txt = fullfile(fp_rep, 'input_set_path.txt');
    end
    fid = fopen(fp_txt, 'w');
    if fid ~= -1
        fprintf(fid, '%s\n', input_file);
        fclose(fid);
    end
catch
end

end
