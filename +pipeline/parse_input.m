function files = parse_input(input_path)
%PARSE_INPUT Resolve EEG .set input files from file/folder/GUI.

files = {};
if isempty(input_path)
    choice = questdlg('Select input type', 'EEG Pipeline', 'File', 'Folder', 'Cancel', 'File');
    if strcmp(choice,'File')
        [fn, fp] = uigetfile('*.set', 'Select .set file');
        if isequal(fn,0)
            error('No file selected');
        end
        files = {fullfile(fp, fn)};
    elseif strcmp(choice,'Folder')
        fp = uigetdir(pwd, 'Select folder containing .set files');
        if isequal(fp,0)
            error('No folder selected');
        end
        d = dir(fullfile(fp, '*.set'));
        if isempty(d)
            error('No .set files found in folder: %s', fp);
        end
        files = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
    else
        error('No input selected');
    end
elseif isfolder(input_path)
    d = dir(fullfile(input_path, '*.set'));
    if isempty(d)
        error('No .set files found in folder: %s', input_path);
    end
    files = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
else
    if ~endsWith(lower(input_path), '.set')
        error('Input file must be a .set file');
    end
    files = {input_path};
end

end
