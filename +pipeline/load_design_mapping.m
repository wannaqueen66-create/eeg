function MapLong = load_design_mapping(design_path)
%LOAD_DESIGN_MAPPING Load and reshape the trial design mapping CSV.
%
% Input CSV format (wide):
%   name, SportFreq, Experience, Order,
%   trial01_scene, trial01_WWR, trial01_Cond, trial01_Complexity,
%   ... trial12_*
%
% Output (long):
%   subject_id, scene_id(1..12), scene_name, WWR, Cond, Complexity,
%   SportFreq, Experience, Order

fp = pipeline.find_design_file(design_path);
if isempty(fp)
    MapLong = table();
    return;
end

T = readtable(fp, 'TextType','string');

% normalize subject id column
if any(strcmpi(T.Properties.VariableNames,'name'))
    subj = T.(T.Properties.VariableNames{strcmpi(T.Properties.VariableNames,'name')});
else
    error('Design file must contain a column named "name" (subject id).');
end

if any(strcmpi(T.Properties.VariableNames,'SportFreq'))
    sport = T.(T.Properties.VariableNames{strcmpi(T.Properties.VariableNames,'SportFreq')});
else
    sport = repmat("", height(T), 1);
end

if any(strcmpi(T.Properties.VariableNames,'Experience'))
    expv = T.(T.Properties.VariableNames{strcmpi(T.Properties.VariableNames,'Experience')});
else
    expv = repmat("", height(T), 1);
end

if any(strcmpi(T.Properties.VariableNames,'Order'))
    orderv = T.(T.Properties.VariableNames{strcmpi(T.Properties.VariableNames,'Order')});
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

        if ~ismember(col_scene, T.Properties.VariableNames)
            error('Missing column %s in design file.', col_scene);
        end

        scene_name = T.(col_scene)(r);
        wwr = NaN; cond = ""; cmp = NaN;
        if ismember(col_WWR,  T.Properties.VariableNames); wwr  = double(T.(col_WWR)(r)); end
        if ismember(col_Cond, T.Properties.VariableNames); cond = string(T.(col_Cond)(r)); end
        if ismember(col_Cmp,  T.Properties.VariableNames); cmp  = double(T.(col_Cmp)(r)); end

        rows = [rows; {sid, s, string(scene_name), wwr, cond, cmp, string(sport(r)), string(expv(r)), double(orderv(r)), string(fp)}]; %#ok<AGROW>
    end
end

MapLong = rows;

end
