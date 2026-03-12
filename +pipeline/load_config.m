function cfg = load_config(config_path)
%LOAD_CONFIG Read JSON config and apply defaults.

cfg = struct();
try
    if exist(config_path, 'file')
        raw = fileread(config_path);
        cfg = jsondecode(raw);
    end
catch
    cfg = struct();
end

% defaults
if ~isfield(cfg,'gray_dur_min'); cfg.gray_dur_min = 3; end
if ~isfield(cfg,'gray_dur_max'); cfg.gray_dur_max = 15; end
if ~isfield(cfg,'quest_dur_min'); cfg.quest_dur_min = 5; end
if ~isfield(cfg,'quest_dur_max'); cfg.quest_dur_max = 120; end
if ~isfield(cfg,'pairing_mode'); cfg.pairing_mode = 'strict'; end
if ~isfield(cfg,'verbose'); cfg.verbose = true; end
if ~isfield(cfg,'log_file'); cfg.log_file = ''; end
if ~isfield(cfg,'output_dir'); cfg.output_dir = ''; end
if ~isfield(cfg,'zip_output'); cfg.zip_output = false; end
if ~isfield(cfg,'global_summary'); cfg.global_summary = false; end
if ~isfield(cfg,'global_summary_path'); cfg.global_summary_path = ''; end
if ~isfield(cfg,'strict_structure'); cfg.strict_structure = true; end
if ~isfield(cfg,'save_log'); cfg.save_log = true; end
if ~isfield(cfg,'figure_visible'); cfg.figure_visible = false; end
if ~isfield(cfg,'batch_summaries'); cfg.batch_summaries = true; end
if ~isfield(cfg,'design_path'); cfg.design_path = ''; end

% Output-layout migration defaults / aliases
if ~isfield(cfg,'timestamp_output_root'); cfg.timestamp_output_root = false; end
if ~isfield(cfg,'timestamp_summary_only'); cfg.timestamp_summary_only = true; end
if ~isfield(cfg,'summary_runs_subdir'); cfg.summary_runs_subdir = 'runs'; end

% Canonical batch-oriented aliases (preferred in main branch docs)
if ~isfield(cfg,'timestamp_batch_runs')
    cfg.timestamp_batch_runs = cfg.timestamp_summary_only;
end
if ~isfield(cfg,'batch_runs_subdir') || isempty(cfg.batch_runs_subdir)
    cfg.batch_runs_subdir = cfg.summary_runs_subdir;
end

end
