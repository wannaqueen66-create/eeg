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

sid = canonical_subject_id_local(subject_id);
map_sid = canonical_subject_id_local(designMap.subject_id);

M = designMap(map_sid==sid, :);
if height(M)==0
    % allow matching when design uses different whitespace/case fallback
    M = designMap(strtrim(lower(map_sid))==strtrim(lower(sid)), :);
end
if height(M)==0
    return;
end

% Join on scene_id (LEFT join: keep Tin rows only)
% NOTE: table variable names must start with a letter (MATLAB identifier rule),
% so do NOT use names like "__orig_row__".
Tin.orig_row_attach_design = (1:height(Tin))';
keep = {'scene_id','scene_name','WWR','Cond','Complexity','SportFreq','Experience','SportFreqGroup','ExperienceGroup','Order', ...
    'Block','Position','Repetition','RepetitionC','SceneID','IPQ_mean','SAM_Valence','Bmean'};
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

function sid = canonical_subject_id_local(x)
sid = strtrim(string(x));
sid = replace(sid, "\", "/");
for i=1:numel(sid)
    parts = split(sid(i), "/");
    sid(i) = parts(end);
end
sid = regexprep(sid, '\\.set$', '', 'ignorecase');
sid = strtrim(sid);
end
