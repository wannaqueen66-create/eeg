function write_chanlocs_snapshot(fp_qc, base, chanlocs)
%WRITE_CHANLOCS_SNAPSHOT Save EEGLAB chanlocs struct for later group topoplots.

if nargin < 3 || isempty(chanlocs)
    return;
end

try
    f = fullfile(fp_qc, sprintf('%s_chanlocs.mat', base));
    save(f, 'chanlocs');
catch
end

end
