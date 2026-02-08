function [fp_sub, fp_csv, fp_fig, fp_qc] = prepare_output(fp, subject_id, cfg)
%PREPARE_OUTPUT Create output directories and write config snapshot.

fp_out = fp;
if isfield(cfg, 'output_dir') && ~isempty(cfg.output_dir)
    fp_out = pipeline.resolve_output_dir(fp, cfg.output_dir);
    if ~exist(fp_out, 'dir')
        mkdir(fp_out);
    end
end

fp_sub = fullfile(fp_out, subject_id);
if ~exist(fp_sub, 'dir')
    mkdir(fp_sub);
end

fp_csv = fullfile(fp_sub, 'csv');
fp_fig = fullfile(fp_sub, 'fig');
fp_qc  = fullfile(fp_sub, 'qc');
if ~exist(fp_csv, 'dir'); mkdir(fp_csv); end
if ~exist(fp_fig, 'dir'); mkdir(fp_fig); end
if ~exist(fp_qc,  'dir'); mkdir(fp_qc); end

pipeline.write_config_snapshot(cfg, fp_sub);

end
