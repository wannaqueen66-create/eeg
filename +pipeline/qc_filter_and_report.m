function out = qc_filter_and_report(fp_out, fp_batch, cfg, AllScene, AllPairs)
%QC_FILTER_AND_REPORT Build QC-based inclusion/exclusion, filter merged tables, and write reports.
%
% Outputs:
%   - qc_exclusion_subjects.csv
%   - qc_scene_valid_counts.csv
%   - all_subjects_scene_level_qc.csv (filtered)
%   - all_subjects_pairs_check_qc.csv (filtered)
%   - qc_filter_report.md (human-readable)

out = struct();

if nargin < 5
    AllPairs = table();
end

if nargin < 4 || isempty(AllScene)
    AllScene = table();
end

if nargin < 3
    cfg = struct();
end

applyQc = true;
try
    if isfield(cfg,'qc_apply'); applyQc = logical(cfg.qc_apply); end
catch
end

hfThr = 0.4;
try
    if isfield(cfg,'qc_hf_threshold'); hfThr = double(cfg.qc_hf_threshold); end
catch
end

badFracThr = 0.3;
try
    if isfield(cfg,'qc_view_bad_frac_threshold'); badFracThr = double(cfg.qc_view_bad_frac_threshold); end
catch
end

rmsMethod = "mad";
try
    if isfield(cfg,'qc_rms_method'); rmsMethod = string(cfg.qc_rms_method); end
catch
end

rmsK = 3.5;
try
    if isfield(cfg,'qc_rms_k'); rmsK = double(cfg.qc_rms_k); end
catch
end

% Read per-subject qc.csv and build view-segment QC map
D = dir(fp_out);
D = D([D.isdir]);
subs = setdiff({D.name},{'.','..','summary'});
subs = sort(subs);

Qsub = table('Size',[0 10], 'VariableTypes', ...
    {'string','double','double','double','double','double','double','double','double','string'}, ...
    'VariableNames', {'subject_id','view_mean_hf','gray_mean_hf','view_mean_rms','gray_mean_rms','view_bad_frac','n_view','n_view_bad','n_gray','notes'});

Qview = table('Size',[0 5], 'VariableTypes', {'string','double','double','double','double'}, ...
    'VariableNames', {'subject_id','scene_id','view_hf_ratio','view_rms_mean','view_near_boundary'});

for k=1:numel(subs)
    sid = string(subs{k});
    f_qc = fullfile(fp_out, char(sid), 'qc', sprintf('%s_qc.csv', sid));
    if ~exist(f_qc,'file')
        continue;
    end

    try
        T = readtable(f_qc, 'TextType','string');
    catch
        try
            T = readtable(f_qc);
        catch
            continue;
        end
    end

    if ~ismember('cond', T.Properties.VariableNames)
        continue;
    end

    cond = string(T.cond);

    % column names (robust)
    hfCol = "hf_ratio_20_40Hz";
    if ~ismember(hfCol, T.Properties.VariableNames)
        % legacy fallback
        cand = T.Properties.VariableNames(contains(T.Properties.VariableNames,'hf_ratio'));
        if ~isempty(cand); hfCol = string(cand{1}); end
    end

    rmsCol = "rms_mean_uV";
    if ~ismember(rmsCol, T.Properties.VariableNames)
        cand = T.Properties.VariableNames(contains(T.Properties.VariableNames,'rms'));
        if ~isempty(cand); rmsCol = string(cand{1}); end
    end

    hf = double(T.(char(hfCol)));
    rms = double(T.(char(rmsCol)));

    isView = cond=="view";
    isGray = cond=="gray";

    n_view = sum(isView);
    n_gray = sum(isGray);

    view_bad = isView & hf > hfThr;
    n_view_bad = sum(view_bad);
    view_bad_frac = nan;
    if n_view>0; view_bad_frac = n_view_bad / n_view; end

    view_mean_hf = mean(hf(isView),'omitnan');
    gray_mean_hf = mean(hf(isGray),'omitnan');
    view_mean_rms = mean(rms(isView),'omitnan');
    gray_mean_rms = mean(rms(isGray),'omitnan');

    notes = "";
    if n_view==0
        notes = notes + "no_view;";
    end

    Qsub = [Qsub; {sid, view_mean_hf, gray_mean_hf, view_mean_rms, gray_mean_rms, view_bad_frac, n_view, n_view_bad, n_gray, notes}]; %#ok<AGROW>

    % Per-scene view QC (for scene-level filtering & counts)
    if ismember('scene_id', T.Properties.VariableNames)
        scene_id = double(T.scene_id);
        nearB = zeros(height(T),1);
        if ismember('near_boundary', T.Properties.VariableNames)
            nearB = double(T.near_boundary);
        end
        keep = isView & ~isnan(scene_id);
        if any(keep)
            Qview = [Qview; table(repmat(sid,sum(keep),1), scene_id(keep), hf(keep), rms(keep), nearB(keep), ...
                'VariableNames', Qview.Properties.VariableNames)]; %#ok<AGROW>
        end
    end
end

% Determine RMS thresholds across subjects (robust)
rmsViewThr = NaN; rmsGrayThr = NaN;
try
    if height(Qsub) > 0
        v = Qsub.view_mean_rms;
        g = Qsub.gray_mean_rms;
        if lower(rmsMethod)=="mad"
            rmsViewThr = median(v,'omitnan') + rmsK * mad(v,1);
            rmsGrayThr = median(g,'omitnan') + rmsK * mad(g,1);
        elseif lower(rmsMethod)=="p95" || lower(rmsMethod)=="percentile"
            rmsViewThr = prctile(v,95);
            rmsGrayThr = prctile(g,95);
        elseif lower(rmsMethod)=="p97_5" || lower(rmsMethod)=="p97.5"
            rmsViewThr = prctile(v,97.5);
            rmsGrayThr = prctile(g,97.5);
        end
    end
catch
end

% Apply exclusion rules (subject-level)
exclude_scene = false(height(Qsub),1);
exclude_recov = false(height(Qsub),1);
reasons_scene = strings(height(Qsub),1);
reasons_recov = strings(height(Qsub),1);

for i=1:height(Qsub)
    r = strings(0,1);
    if ~isnan(Qsub.view_mean_hf(i)) && Qsub.view_mean_hf(i) > hfThr; r(end+1) = "view_mean_hf"; end
    if ~isnan(Qsub.view_bad_frac(i)) && Qsub.view_bad_frac(i) >= badFracThr; r(end+1) = "view_bad_frac"; end
    if ~isnan(rmsViewThr) && ~isnan(Qsub.view_mean_rms(i)) && Qsub.view_mean_rms(i) > rmsViewThr; r(end+1) = "view_mean_rms"; end

    exclude_scene(i) = applyQc && ~isempty(r);
    reasons_scene(i) = strjoin(r, ";");

    rr = r;
    if ~isnan(Qsub.gray_mean_hf(i)) && Qsub.gray_mean_hf(i) > hfThr; rr(end+1) = "gray_mean_hf"; end
    if ~isnan(rmsGrayThr) && ~isnan(Qsub.gray_mean_rms(i)) && Qsub.gray_mean_rms(i) > rmsGrayThr; rr(end+1) = "gray_mean_rms"; end

    exclude_recov(i) = applyQc && ~isempty(rr);
    reasons_recov(i) = strjoin(rr, ";");
end

Qsub.exclude_scenelevel = exclude_scene;
Qsub.exclude_recovery = exclude_recov;
Qsub.reasons_scenelevel = reasons_scene;
Qsub.reasons_recovery = reasons_recov;
Qsub.rms_view_threshold = repmat(rmsViewThr, height(Qsub), 1);
Qsub.rms_gray_threshold = repmat(rmsGrayThr, height(Qsub), 1);

fp_tbl_qc = fp_batch;
fp_rep = fp_batch;
fp_aud = fp_batch;
try
    if exist(fullfile(fp_batch, 'qc'),'dir')
        fp_tbl_qc = fullfile(fp_batch, 'qc');
    elseif exist('pipeline.get_table_dir','file')==2
        fp_tbl_qc = pipeline.get_table_dir(fp_batch, cfg, 'merged_qc');
    end
    if exist(fullfile(fp_batch, 'reports'),'dir')
        fp_rep = fullfile(fp_batch, 'reports');
    elseif exist('pipeline.get_report_dir','file')==2
        fp_rep = pipeline.get_report_dir(fp_batch, cfg);
    end
    if exist(fullfile(fp_batch, 'audit'),'dir')
        fp_aud = fullfile(fp_batch, 'audit');
    elseif exist('pipeline.get_audit_dir','file')==2
        fp_aud = pipeline.get_audit_dir(fp_batch, cfg);
    end
catch
end

out.qc_exclusion_subjects = fullfile(fp_tbl_qc, 'qc_exclusion_subjects.csv');
try
    writetable(Qsub, out.qc_exclusion_subjects);
    fprintf('[QC] Wrote subject exclusion table: %s\n', out.qc_exclusion_subjects);
catch
end

% Filter scene-level rows using:
%   - subject-level exclusion (exclude_scenelevel)
%   - per-scene view HF_ratio (segment-level): hf > hfThr
AllScene_qc = AllScene;
if ~isempty(AllScene_qc) && ismember('subject_id', AllScene_qc.Properties.VariableNames)
    % join subject exclusion
    [tf, loc] = ismember(string(AllScene_qc.subject_id), string(Qsub.subject_id));
    exSub = false(height(AllScene_qc),1);
    exSub(tf) = logical(Qsub.exclude_scenelevel(loc(tf)));

    % join per-scene view hf/rms
    if ~isempty(Qview) && ismember('scene_id', AllScene_qc.Properties.VariableNames)
        key1 = string(AllScene_qc.subject_id) + "__" + string(AllScene_qc.scene_id);
        key2 = string(Qview.subject_id) + "__" + string(Qview.scene_id);
        [tf2, loc2] = ismember(key1, key2);
        view_hf = nan(height(AllScene_qc),1);
        view_rms = nan(height(AllScene_qc),1);
        if any(tf2)
            view_hf(tf2) = Qview.view_hf_ratio(loc2(tf2));
            view_rms(tf2) = Qview.view_rms_mean(loc2(tf2));
        end
        AllScene_qc.view_hf_ratio = view_hf;
        AllScene_qc.view_rms_mean = view_rms;
    end

    badSeg = false(height(AllScene_qc),1);
    if ismember('view_hf_ratio', AllScene_qc.Properties.VariableNames)
        badSeg = AllScene_qc.view_hf_ratio > hfThr;
    end

    keep = ~exSub & ~badSeg;
    out.n_scene_rows_before = height(AllScene_qc);
    out.n_scene_rows_after = sum(keep);
    AllScene_qc = AllScene_qc(keep,:);
end

% Filter pair-level recovery rows using:
%   - subject-level exclusion for recovery
AllPairs_qc = AllPairs;
if ~isempty(AllPairs_qc) && ismember('subject_id', AllPairs_qc.Properties.VariableNames)
    [tf, loc] = ismember(string(AllPairs_qc.subject_id), string(Qsub.subject_id));
    exSub = false(height(AllPairs_qc),1);
    exSub(tf) = logical(Qsub.exclude_recovery(loc(tf)));

    % also drop pairs whose view scene had hf>thr (segment-level)
    badSeg = false(height(AllPairs_qc),1);
    if ~isempty(Qview) && ismember('scene_id', AllPairs_qc.Properties.VariableNames)
        key1 = string(AllPairs_qc.subject_id) + "__" + string(AllPairs_qc.scene_id);
        key2 = string(Qview.subject_id) + "__" + string(Qview.scene_id);
        [tf2, loc2] = ismember(key1, key2);
        if any(tf2)
            badSeg(tf2) = Qview.view_hf_ratio(loc2(tf2)) > hfThr;
        end
    end

    keep = ~exSub & ~badSeg;
    out.n_pairs_rows_before = height(AllPairs_qc);
    out.n_pairs_rows_after = sum(keep);
    AllPairs_qc = AllPairs_qc(keep,:);
end

% Write filtered merged tables
if ~isempty(AllScene_qc)
    out.all_subjects_scene_level_qc = fullfile(fp_tbl_qc, 'all_subjects_scene_level_qc.csv');
    try; writetable(AllScene_qc, out.all_subjects_scene_level_qc); catch; end
end
if ~isempty(AllPairs_qc)
    out.all_subjects_pairs_check_qc = fullfile(fp_tbl_qc, 'all_subjects_pairs_check_qc.csv');
    try; writetable(AllPairs_qc, out.all_subjects_pairs_check_qc); catch; end
end

% Scene valid counts + excluded subject names
% Important semantic rule:
% - scene_id is a presentation/trial index (Block+Position), not guaranteed to be
%   a stable scene identity across counterbalanced orders.
% - If the same scene_id maps to multiple scene_name values across subjects,
%   aggregate by scene_name instead and mark the output accordingly.
try
    if ~isempty(AllScene) && ismember('scene_id', AllScene.Properties.VariableNames) && ismember('subject_id', AllScene.Properties.VariableNames)
        sc = AllScene;
        sc.subject_id = string(sc.subject_id);
        exRow = strings(height(sc),1);

        [tf, loc] = ismember(string(sc.subject_id), string(Qsub.subject_id));
        exSub = false(height(sc),1);
        exSub(tf) = logical(Qsub.exclude_scenelevel(loc(tf)));
        exRow(exSub) = "subject_excluded";

        if ~isempty(Qview)
            key1 = string(sc.subject_id) + "__" + string(sc.scene_id);
            key2 = string(Qview.subject_id) + "__" + string(Qview.scene_id);
            [tf2, loc2] = ismember(key1, key2);
            hf = nan(height(sc),1);
            if any(tf2); hf(tf2) = Qview.view_hf_ratio(loc2(tf2)); end
            exHf = hf > hfThr;
            exRow(exHf & exRow=="") = "view_hf_excluded";
        end

        sc.ex_reason = exRow;

        useSceneIdentity = false;
        if ismember('scene_name', sc.Properties.VariableNames)
            tmp = sc(:, {'scene_id','scene_name'});
            tmp.scene_name = strtrim(string(tmp.scene_name));
            tmp = tmp(strlength(tmp.scene_name)>0 & ~isnan(tmp.scene_id), :);
            if ~isempty(tmp)
                Gtmp = groupsummary(tmp, 'scene_id', 'numel', 'scene_name');
                if ismember('GroupCount', Gtmp.Properties.VariableNames)
                    useSceneIdentity = any(double(Gtmp.GroupCount) > 1);
                end
            end
        end

        rows = {};
        if useSceneIdentity && ismember('scene_name', sc.Properties.VariableNames)
            ids = unique(strtrim(string(sc.scene_name)), 'stable');
            ids = ids(strlength(ids)>0);
            for ii = 1:numel(ids)
                sidName = ids(ii);
                idx = strtrim(string(sc.scene_name)) == sidName;
                n_total = sum(idx);
                n_valid = sum(idx & sc.ex_reason=="");
                excl = sc.subject_id(idx & sc.ex_reason~="");
                excl = unique(excl,'stable');
                rows(end+1,:) = {sidName, n_total, n_valid, strjoin(excl, ';')}; %#ok<AGROW>
            end
            Tcnt = cell2table(rows, 'VariableNames', {'scene_identity','n_total','n_valid','excluded_subjects'});
        else
            sceneU = unique(sc.scene_id);
            sceneU = sort(sceneU(:));
            for si = sceneU'
                idx = sc.scene_id==si;
                n_total = sum(idx);
                n_valid = sum(idx & sc.ex_reason=="");
                excl = sc.subject_id(idx & sc.ex_reason~="");
                excl = unique(excl,'stable');
                rows(end+1,:) = {double(si), n_total, n_valid, strjoin(excl, ';')}; %#ok<AGROW>
            end
            Tcnt = cell2table(rows, 'VariableNames', {'scene_id','n_total','n_valid','excluded_subjects'});
            if ismember('scene_name', sc.Properties.VariableNames)
                labs = strings(height(Tcnt),1);
                for i=1:height(Tcnt)
                    si = Tcnt.scene_id(i);
                    idx = sc.scene_id==si;
                    sn = "";
                    try
                        sn = string(sc.scene_name(find(idx,1,'first')));
                    catch
                    end
                    labs(i) = sn;
                end
                Tcnt.scene_name = labs;
            end
        end

        out.qc_scene_valid_counts = fullfile(fp_tbl_qc, 'qc_scene_valid_counts.csv');
        writetable(Tcnt, out.qc_scene_valid_counts);
        fprintf('[QC] Wrote scene valid counts: %s\n', out.qc_scene_valid_counts);
    end
catch
end

% Human-readable report
try
    fp_md = fullfile(fp_rep, 'qc_filter_report.md');
    fid = fopen(fp_md,'w');
    if fid~=-1
        fprintf(fid, '# QC Filter Report (Auto-generated)\n\n');
        fprintf(fid, '- Generated: %s\n', datestr(now,31));
        fprintf(fid, '- HF threshold (20-40 / 1-40): %.3f\n', hfThr);
        fprintf(fid, '- View bad fraction threshold: %.3f\n', badFracThr);
        fprintf(fid, '- RMS method: %s\n', rmsMethod);
        fprintf(fid, '- RMS k: %.2f\n', rmsK);
        fprintf(fid, '- RMS view threshold (subject mean): %.3f µV\n', rmsViewThr);
        fprintf(fid, '- RMS gray threshold (subject mean): %.3f µV\n\n', rmsGrayThr);

        if height(Qsub)>0
            fprintf(fid, '## Subject exclusions\n\n');
            fprintf(fid, '- n_subjects: %d\n', height(Qsub));
            fprintf(fid, '- excluded (scene-level): %d\n', sum(Qsub.exclude_scenelevel));
            fprintf(fid, '- excluded (recovery): %d\n\n', sum(Qsub.exclude_recovery));
        end

        % Scene-level valid counts (after exclusions)
        try
            if isfield(out,'qc_scene_valid_counts') && exist(out.qc_scene_valid_counts,'file')
                Tcnt = readtable(out.qc_scene_valid_counts, 'TextType','string');
                fprintf(fid, '## Scene valid counts (after QC)\n\n');
                if ismember('scene_identity', Tcnt.Properties.VariableNames)
                    fprintf(fid, '| scene_identity | valid / total | excluded_subjects |\n');
                    fprintf(fid, '|---|---:|---|\n');
                    for i=1:height(Tcnt)
                        sidName = string(Tcnt.scene_identity(i));
                        excl = "";
                        if ismember('excluded_subjects', Tcnt.Properties.VariableNames)
                            excl = string(Tcnt.excluded_subjects(i));
                        end
                        fprintf(fid, '| %s | %d / %d | %s |\n', sidName, Tcnt.n_valid(i), Tcnt.n_total(i), excl);
                    end
                else
                    fprintf(fid, '| trial_index | scene_name | valid / total | excluded_subjects |\n');
                    fprintf(fid, '|---:|---|---:|---|\n');
                    for i=1:height(Tcnt)
                        sid = Tcnt.scene_id(i);
                        sname = "";
                        if ismember('scene_name', Tcnt.Properties.VariableNames)
                            sname = string(Tcnt.scene_name(i));
                        end
                        excl = "";
                        if ismember('excluded_subjects', Tcnt.Properties.VariableNames)
                            excl = string(Tcnt.excluded_subjects(i));
                        end
                        fprintf(fid, '| %d | %s | %d / %d | %s |\n', sid, sname, Tcnt.n_valid(i), Tcnt.n_total(i), excl);
                    end
                end
                fprintf(fid, '\n');
            end
        catch
        end

        fclose(fid);
    end
    out.qc_filter_report_md = fp_md;
catch
end

out.AllScene_qc = AllScene_qc;
out.AllPairs_qc = AllPairs_qc;
out.Qsub = Qsub;
out.Qview = Qview;

% Write machine-readable QC audit trail JSON (for supplementary materials)
try
    doJson = true;
    if isfield(cfg,'qc_write_json'); doJson = logical(cfg.qc_write_json); end
    if doJson
        qc = struct();
        qc.generated = datestr(now,31);
        qc.git_commit = '';
        try
            [st,outc] = system('git rev-parse HEAD');
            if st==0; qc.git_commit = strtrim(outc); end
        catch
        end

        qc.params = struct();
        qc.params.qc_apply = applyQc;
        qc.params.qc_hf_threshold = hfThr;
        qc.params.qc_view_bad_frac_threshold = badFracThr;
        qc.params.qc_rms_method = char(rmsMethod);
        qc.params.qc_rms_k = rmsK;

        qc.derived = struct();
        qc.derived.rms_view_threshold_subject_mean = rmsViewThr;
        qc.derived.rms_gray_threshold_subject_mean = rmsGrayThr;

        qc.summary = struct();
        qc.summary.n_subjects = height(Qsub);
        qc.summary.n_excluded_scenelevel = sum(Qsub.exclude_scenelevel);
        qc.summary.n_excluded_recovery = sum(Qsub.exclude_recovery);

        qc.excluded = struct();
        qc.excluded.scenelevel_subjects = cellstr(string(Qsub.subject_id(Qsub.exclude_scenelevel)));
        qc.excluded.recovery_subjects = cellstr(string(Qsub.subject_id(Qsub.exclude_recovery)));

        % Per-scene valid counts (as struct array)
        qc.scene_counts = [];
        try
            if isfield(out,'qc_scene_valid_counts') && exist(out.qc_scene_valid_counts,'file')
                Tcnt = readtable(out.qc_scene_valid_counts, 'TextType','string');
                arr = repmat(struct('scene_id',[],'scene_name','','n_total',[],'n_valid',[],'excluded_subjects',{{}}), height(Tcnt), 1);
                for i=1:height(Tcnt)
                    arr(i).scene_id = double(Tcnt.scene_id(i));
                    if ismember('scene_name', Tcnt.Properties.VariableNames)
                        arr(i).scene_name = char(string(Tcnt.scene_name(i)));
                    end
                    arr(i).n_total = double(Tcnt.n_total(i));
                    arr(i).n_valid = double(Tcnt.n_valid(i));
                    excl = '';
                    if ismember('excluded_subjects', Tcnt.Properties.VariableNames)
                        excl = char(string(Tcnt.excluded_subjects(i)));
                    end
                    if strlength(string(excl))>0
                        arr(i).excluded_subjects = strsplit(excl, ';');
                    else
                        arr(i).excluded_subjects = {};
                    end
                end
                qc.scene_counts = arr;
            end
        catch
        end

        fp_json = fullfile(fp_aud, 'qc_report.json');
        fid = fopen(fp_json,'w');
        if fid~=-1
            fwrite(fid, jsonencode(qc, 'PrettyPrint', true), 'char');
            fclose(fid);
            out.qc_report_json = fp_json;
            fprintf('[QC] Wrote QC audit JSON: %s\n', fp_json);
        end
    end
catch
end

end
