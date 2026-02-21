function fp = find_design_file(design_path)
%FIND_DESIGN_FILE Resolve a design CSV file path.
% Accepts either a direct CSV path or a folder containing a CSV.

fp = '';
if nargin < 1 || isempty(design_path)
    return;
end

if isfolder(design_path)
    % Prefer common filenames
    candidates = { ...
        fullfile(design_path,'scene_mapping.csv'), ...
        fullfile(design_path,'design.csv'), ...
        fullfile(design_path,'mapping.csv') ...
    };
    for i = 1:numel(candidates)
        if exist(candidates{i},'file')
            fp = candidates{i};
            return;
        end
    end
    d = dir(fullfile(design_path,'*.csv'));
    if ~isempty(d)
        fp = fullfile(d(1).folder, d(1).name);
    end
else
    if exist(design_path,'file')
        fp = design_path;
    end
end

end
