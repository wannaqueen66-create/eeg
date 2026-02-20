function s = sanitize_filename(s)
%SANITIZE_FILENAME Make a string safe for use as a filename on Windows.

s = string(s);
% Replace invalid characters
s = regexprep(s, '[<>:"/\\|?*]', '_');
% Trim trailing spaces/dots (Windows)
s = regexprep(s, '[\s\.]+$', '');
if strlength(s) == 0
    s = "untitled";
end
s = char(s);

end
