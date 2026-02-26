function fp_out = get_output_root(fp_in, cfg)
%GET_OUTPUT_ROOT Return the root output directory for a given input folder.
%
% Default output root:
%   <fp_in>/bandpower_outputs
%
% Optional timestamped output root (archive runs without overwrite):
%   If cfg.timestamp_output_root==true, output root becomes:
%     <fp_in>/bandpower_outputs_YYYYMMDD_HHMMSS
%   and we write/update a pointer file:
%     <fp_in>/bandpower_outputs_latest.txt
%   so you can find the latest run easily.
%
% If cfg.output_dir is set, it takes precedence.

% 1) Explicit output_dir wins
if isfield(cfg, 'output_dir') && ~isempty(cfg.output_dir)
    fp_out = pipeline.resolve_output_dir(fp_in, cfg.output_dir);
else
    base = 'bandpower_outputs';

    useTs = false;
    try
        if isfield(cfg,'timestamp_output_root')
            useTs = logical(cfg.timestamp_output_root);
        end
    catch
        useTs = false;
    end

    if useTs
        % Keep a stable timestamp within this MATLAB session
        persistent RUN_TS;
        if isempty(RUN_TS)
            RUN_TS = datestr(now, 'yyyymmdd_HHMMSS');
        end
        outName = sprintf('%s_%s', base, RUN_TS);
        fp_out = fullfile(fp_in, outName);

        % Write/update latest pointer
        try
            fp_ptr = fullfile(fp_in, 'bandpower_outputs_latest.txt');
            fid = fopen(fp_ptr, 'w');
            if fid > 0
                fprintf(fid, '%s\n', fp_out);
                fclose(fid);
            end
        catch
        end
    else
        fp_out = fullfile(fp_in, base);
    end
end

if ~exist(fp_out, 'dir')
    mkdir(fp_out);
end

end
