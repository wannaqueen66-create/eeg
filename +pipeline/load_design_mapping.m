function MapLong = load_design_mapping(design_path)
%LOAD_DESIGN_MAPPING Load and reshape the trial design mapping CSV.
%
% Input CSV format supports BOTH:
%
% (A) Wide design (legacy):
%   name, SportFreq, Experience, Order,
%   trial01_scene, trial01_WWR, trial01_Cond, trial01_Complexity, ... trial12_*
%
% (B) Long design (recommended, supports per-scene ratings):
%   SubjectID, Order, Block, Position, WWR, Condition, Complexity,
%   SportFreq, Experience, (optional) ExperienceGroup, SportFreqGroup,
%   (optional) S1.., Bmean, SceneID, ...
%
% Output (long):
%   subject_id, scene_id(1..12), scene_name, WWR, Cond, Complexity,
%   SportFreq, Experience, Order, ... plus any extra columns from design.

fp = pipeline.find_design_file(design_path);
if isempty(fp)
    MapLong = table();
    return;
end

T = readtable(fp, 'TextType','string');
vars = T.Properties.VariableNames;

% --- Robust header cleanup (handles UTF-8 BOM like '\ufeffSubjectID') ---
vars_clean = regexprep(vars, '^\x{FEFF}', '');
if any(~strcmp(vars, vars_clean))
    T.Properties.VariableNames = vars_clean;
    vars = vars_clean;
end

% Detect format
isWide = ismember('trial01_scene', vars) || any(startsWith(lower(vars), 'trial01_'));
% Long format only requires SubjectID + Block + Position; SceneID is optional
isLong = ismember('SubjectID', vars) && ismember('Block', vars) && ismember('Position', vars);

if isLong && ~isWide
    % ---------- Long format ----------
    % Required columns
    if ~ismember('SubjectID', vars); error('Long design file missing SubjectID'); end
    if ~ismember('Block', vars);     error('Long design file missing Block'); end
    if ~ismember('Position', vars);  error('Long design file missing Position'); end

    subj = strtrim(T.SubjectID);
    block_id = double(T.Block);
    cycle_in_block = double(T.Position);
    scene_id = (block_id-1).*6 + cycle_in_block;

    % Scene identity
    scene_name = "";
    if ismember('SceneID', vars)
        scene_name = string(T.SceneID);
    end

    % Standard factors
    wwr = NaN(height(T),1);
    if ismember('WWR', vars); wwr = double(T.WWR); end

    cond = "";
    if ismember('Condition', vars); cond = string(T.Condition); end

    cmp = NaN(height(T),1);
    if ismember('Complexity', vars); cmp = double(T.Complexity); end

    sport = "";
    % Prefer explicit High/Low grouping column when present.
    if ismember('SportFreqGroup', vars)
        sport = string(T.SportFreqGroup);
    elseif ismember('SportFreq', vars)
        sport = string(T.SportFreq);
    end

    expv = "";
    % Prefer explicit High/Low grouping column when present.
    if ismember('ExperienceGroup', vars)
        expv = string(T.ExperienceGroup);
    elseif ismember('Experience', vars)
        expv = string(T.Experience);
    end

    % If both score and group columns exist but selected value is numeric/empty,
    % fallback to the explicit group column.
    if ismember('SportFreqGroup', vars)
        tmp = strtrim(lower(sport));
        bad = strlength(tmp)==0 | (~ismember(tmp,["high","low","高","低","1","0","h","l"]));
        if any(bad)
            sport(bad) = string(T.SportFreqGroup(bad));
        end
    end
    if ismember('ExperienceGroup', vars)
        tmp = strtrim(lower(expv));
        bad = strlength(tmp)==0 | (~ismember(tmp,["high","low","高","低","1","0","h","l"]));
        if any(bad)
            expv(bad) = string(T.ExperienceGroup(bad));
        end
    end

    orderv = nan(height(T),1);
    if ismember('Order', vars); orderv = double(T.Order); end

    sport_grp = repmat("", height(T), 1);
    exp_grp = repmat("", height(T), 1);
    if ismember('SportFreqGroup', vars); sport_grp = string(T.SportFreqGroup); end
    if ismember('ExperienceGroup', vars); exp_grp = string(T.ExperienceGroup); end

    rows = table(subj, scene_id, scene_name, wwr, cond, cmp, sport, expv, sport_grp, exp_grp, orderv, repmat(string(fp),height(T),1), ...
        'VariableNames', {'subject_id','scene_id','scene_name','WWR','Cond','Complexity','SportFreq','Experience','SportFreqGroup','ExperienceGroup','Order','design_file'});

    % Attach extra columns (ratings etc.)
    extraVars = setdiff(vars, {'SubjectID','Block','Position','WWR','Condition','Complexity','SportFreq','Experience','ExperienceGroup','SportFreqGroup','Order','SceneID'});
    for i = 1:numel(extraVars)
        v = extraVars{i};
        % Avoid name collisions
        if ismember(v, rows.Properties.VariableNames)
            rows.("design_"+v) = T.(v);
        else
            rows.(v) = T.(v);
        end
    end

    MapLong = rows;
    return;
end

% ---------- Wide format (legacy) ----------
% normalize subject id column
if any(strcmpi(vars,'name'))
    subj = T.(vars{strcmpi(vars,'name')});
else
    error('Wide design file must contain a column named "name" (subject id).');
end

if any(strcmpi(vars,'SportFreq'))
    sport = T.(vars{strcmpi(vars,'SportFreq')});
else
    sport = repmat("", height(T), 1);
end

if any(strcmpi(vars,'Experience'))
    expv = T.(vars{strcmpi(vars,'Experience')});
else
    expv = repmat("", height(T), 1);
end

if any(strcmpi(vars,'Order'))
    orderv = T.(vars{strcmpi(vars,'Order')});
else
    orderv = nan(height(T), 1);
end

rows = table('Size',[0 10], ...
    'VariableTypes', {'string','double','string','double','string','double','string','string','double','string'}, ...
    'VariableNames', {'subject_id','scene_id','scene_name','WWR','Cond','Complexity','SportFreq','Experience','Order','design_file'});

for r = 1:height(T)
    sid = strtrim(subj(r));
    for s = 1:12
        col_scene = sprintf('trial%02d_scene', s);
        col_WWR   = sprintf('trial%02d_WWR', s);
        col_Cond  = sprintf('trial%02d_Cond', s);
        col_Cmp   = sprintf('trial%02d_Complexity', s);

        if ~ismember(col_scene, vars)
            error('Missing column %s in design file.', col_scene);
        end

        scene_name = T.(col_scene)(r);
        wwr = NaN; cond = ""; cmp = NaN;
        if ismember(col_WWR,  vars); wwr  = double(T.(col_WWR)(r)); end
        if ismember(col_Cond, vars); cond = string(T.(col_Cond)(r)); end
        if ismember(col_Cmp,  vars); cmp  = double(T.(col_Cmp)(r)); end

        rows = [rows; {sid, s, string(scene_name), wwr, cond, cmp, string(sport(r)), string(expv(r)), double(orderv(r)), string(fp)}]; %#ok<AGROW>
    end
end

MapLong = rows;

end
