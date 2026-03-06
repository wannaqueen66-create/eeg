function [Tout, meta] = merge_eye_scene_features(Tscene, eye_summary_path, cfg)
%MERGE_EYE_SCENE_FEATURES Merge scene-level eye-tracking features into EEG scene table.
%
% Inputs
%   Tscene            EEG merged/all-subjects scene-level table
%   eye_summary_path  path to eye scene-level CSV
%   cfg               config struct (optional)
%
% Expected eye CSV keys (minimum):
%   - subject_id
%   - scene_id
%
% Output
%   Tout              merged table (left join on EEG rows)
%   meta              struct with merge diagnostics
%
% Notes
%   - Eye feature columns are auto-prefixed with `eye_` when needed to avoid
%     collisions with existing EEG columns.
%   - Join keys default to {'subject_id','scene_id'} unless overridden by
%     cfg.eye_merge_keys.

Tout = Tscene;
meta = struct('ok',false,'path',string(eye_summary_path), ...
    'n_eeg_rows',height(Tscene),'n_eye_rows',0,'n_merged_rows',height(Tscene), ...
    'keys',{{'subject_id','scene_id'}},'matched_rows',0,'missing_eye_rows',height(Tscene), ...
    'renamed_columns',strings(0,1));

if nargin < 2 || isempty(eye_summary_path) || ~isfile(eye_summary_path)
    return;
end
if nargin < 3
    cfg = struct();
end
if isempty(Tscene) || ~istable(Tscene)
    return;
end

keys = {'subject_id','scene_id'};
try
    if isfield(cfg,'eye_merge_keys') && ~isempty(cfg.eye_merge_keys)
        keys = cellstr(string(cfg.eye_merge_keys));
    end
catch
end
meta.keys = keys;

Eye = readtable(eye_summary_path, 'TextType','string');
meta.n_eye_rows = height(Eye);
if isempty(Eye)
    return;
end

for i = 1:numel(keys)
    if ~ismember(keys{i}, Tscene.Properties.VariableNames)
        warning('merge_eye_scene_features: EEG table missing key: %s', keys{i});
        return;
    end
    if ~ismember(keys{i}, Eye.Properties.VariableNames)
        warning('merge_eye_scene_features: eye table missing key: %s', keys{i});
        return;
    end
end

% Canonicalize key columns
for i = 1:numel(keys)
    k = keys{i};
    if iscellstr(Tscene.(k)) || isstring(Tscene.(k)) || ischar(Tscene.(k))
        Tout.(k) = strtrim(string(Tscene.(k)));
    end
    if iscellstr(Eye.(k)) || isstring(Eye.(k)) || ischar(Eye.(k))
        Eye.(k) = strtrim(string(Eye.(k)));
    end
end

% Prefix conflicting/non-key columns with eye_
eyeVars = Eye.Properties.VariableNames;
renamed = strings(0,1);
for i = 1:numel(eyeVars)
    v = eyeVars{i};
    if ismember(v, keys)
        continue;
    end
    newv = v;
    if ismember(v, Tout.Properties.VariableNames) || ~startsWith(string(v), "eye_")
        newv = char("eye_" + string(v));
    end
    if ~strcmp(v, newv)
        % ensure uniqueness
        basev = string(newv);
        n = 1;
        while ismember(char(newv), Eye.Properties.VariableNames) || ismember(char(newv), Tout.Properties.VariableNames)
            if strcmp(char(newv), v)
                newv = char(basev + "_" + string(n));
            else
                if ~ismember(char(newv), Eye.Properties.VariableNames) && ~ismember(char(newv), Tout.Properties.VariableNames)
                    break;
                end
                n = n + 1;
                newv = char(basev + "_" + string(n));
            end
        end
        Eye.Properties.VariableNames{strcmp(Eye.Properties.VariableNames, v)} = newv;
        renamed(end+1,1) = string(v) + "->" + string(newv); %#ok<AGROW>
    end
end
meta.renamed_columns = renamed;

% Remove duplicate eye rows by keys (keep first)
try
    [~, ia] = unique(Eye(:, keys), 'rows', 'stable');
    Eye = Eye(ia, :);
catch
end

Tout = outerjoin(Tout, Eye, 'Keys', keys, 'MergeKeys', true, 'Type', 'left');
meta.n_merged_rows = height(Tout);

% Estimate matched rows by presence of first non-key eye variable
nonKeys = setdiff(Eye.Properties.VariableNames, keys, 'stable');
if ~isempty(nonKeys)
    probe = nonKeys{1};
    try
        x = Tout.(probe);
        if isstring(x)
            matched = sum(strlength(x) > 0);
        elseif isnumeric(x)
            matched = sum(~isnan(x));
        elseif iscell(x)
            matched = sum(~cellfun(@isempty, x));
        else
            matched = sum(~ismissing(x));
        end
        meta.matched_rows = matched;
        meta.missing_eye_rows = height(Tout) - matched;
    catch
    end
end

meta.ok = true;
end
