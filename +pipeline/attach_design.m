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

% Join on scene_id
J = outerjoin(Tin, M(:,{'scene_id','scene_name','WWR','Cond','Complexity','SportFreq','Experience','Order'}), ...
    'Keys', 'scene_id', 'MergeKeys', true);

% Restore original row order
if ismember('scene_id', J.Properties.VariableNames)
    % outerjoin may reorder; keep stable by sorting on scene_id then original index if present
end

Tout = J;

end
