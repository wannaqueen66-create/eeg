function Tout = attach_design(Tin, subject_id, designMap)
%ATTACH_DESIGN Attach scene mapping + subject grouping columns to a table.
%
% Tin must contain a numeric scene_id column.
% subject_id: base filename (no extension)
% designMap: output of pipeline.load_design_mapping

Tout = Tin;
if isempty(designMap) || ~istable(designMap) || height(designMap)==0
    return;
end
if ~ismember('scene_id', Tin.Properties.VariableNames)
    return;
end

sid = string(subject_id);

M = designMap(designMap.subject_id==sid, :);
if height(M)==0
    % allow matching when design uses different whitespace
    M = designMap(strtrim(designMap.subject_id)==strtrim(sid), :);
end
if height(M)==0
    return;
end

% Join on scene_id (LEFT join: keep Tin rows only)
% NOTE: table variable names must start with a letter (MATLAB identifier rule),
% so do NOT use names like "__orig_row__".
Tin.orig_row_attach_design = (1:height(Tin))';
keep = {'scene_id','scene_name','WWR','Cond','Complexity','SportFreq','Experience','SportFreqGroup','ExperienceGroup','Order'};
keep = keep(ismember(keep, M.Properties.VariableNames));
J = outerjoin(Tin, M(:,keep), ...
    'Keys', 'scene_id', 'MergeKeys', true, 'Type','left');

% Restore original row order
if ismember('orig_row_attach_design', J.Properties.VariableNames)
    J = sortrows(J, 'orig_row_attach_design');
    J.orig_row_attach_design = [];
end

Tout = J;

end
