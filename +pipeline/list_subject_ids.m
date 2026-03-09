function subs = list_subject_ids(fp_in, cfg)
%LIST_SUBJECT_IDS Return subject ids from the active subject output root.
%
% Stage-2 refactor behavior:
% - Prefer the new subject directory root when available
% - Fallback to historical output-root subject folders for compatibility

subs = strings(0,1);

% New layout
try
    fp_run = pipeline.get_run_dir(fp_in, cfg);
    fp_subjects = fullfile(fp_run, 'subjects');
    if exist(fp_subjects,'dir')
        D = dir(fp_subjects);
        D = D([D.isdir]);
        names = setdiff({D.name}, {'.','..'});
        if ~isempty(names)
            subs = string(sort(names));
            return;
        end
    end
catch
end

% Legacy fallback
try
    fp_out = pipeline.get_output_root(fp_in, cfg);
    D = dir(fp_out);
    D = D([D.isdir]);
    names = setdiff({D.name}, {'.','..','summary','runs'});
    subs = string(sort(names));
catch
    subs = strings(0,1);
end
end
