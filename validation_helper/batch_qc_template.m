function qc = batch_qc_template(conditioned_dir, metrics_dir, species)
%BATCH_QC_TEMPLATE  Reference pattern for wiring QC gates into batch processing.
%
%   qc = batch_qc_template(conditioned_dir, metrics_dir, species)
%
%   This is a TEMPLATE, not a finished pipeline — adapt the two "STAGE" blocks
%   to call your actual conditioning / feature-extraction entry points.  The
%   point it demonstrates: every stage validates its OWN output as it is
%   produced, records accumulate, and qc_summarize emits one report at the end.
%
%   The validators are cheap and side-effect-free, so running them on every
%   recording costs almost nothing and turns silent batch failures into a
%   reviewable list.
%
%   Inputs
%     conditioned_dir  folder of conditioned .mat files (one per recording).
%     metrics_dir      folder where feature/metrics .mat files are written.
%     species          "rat" | "rabbit" | "pig" | "human" (for DF-band hints).
%
%   Output
%     qc               struct with .records (all qc_check records),
%                      .summary and .detail tables from qc_summarize.

    if nargin < 3, species = ""; end

    files = dir(fullfile(conditioned_dir, '*.mat'));
    recs  = repmat(qc_check("","","" ,"PASS",NaN,""), 0, 1);

    for i = 1:numel(files)
        name = files(i).name;
        fprintf('[%d/%d] %s\n', i, numel(files), name);

        % ===== Load + schema check at the boundary ============================
        % "warn" (not "error") so a malformed file is recorded and the batch
        % continues. The schema records go into both copies like any other.
        [cmos_all_data, sr] = load_cmos(fullfile(conditioned_dir, name), ...
                                        "conditioned", "warn");
        cmos_all_data.qc = qc_attach(existing_qc(cmos_all_data), sr);
        recs = [recs; sr]; %#ok<AGROW>
        if any_fail(sr, name, "schema")
            fprintf('  schema FAILED — skipping recording.\n');
            continue;
        end

        % ===== STAGE: conditioning QC =========================================
        % (Already-conditioned file on disk; validate what was produced.)
        % Produce the records ONCE, then write both copies:
        %   (1) embed in cmos_all_data.qc  — provenance, read by the analysis app
        %   (2) append to the batch aggregate — cheap source for qc_summarize
        r = validate_conditioned(cmos_all_data, name);
        cmos_all_data.qc = qc_attach(existing_qc(cmos_all_data), r);
        recs = [recs; r]; %#ok<AGROW>

        % Gate: skip downstream work on a FAILed recording so you don't waste
        % compute extracting features from untrustworthy data.
        if any_fail(r, name, "conditioning")
            fprintf('  conditioning FAILED — skipping feature extraction.\n');
            % Persist the embedded QC even on failure so the file records why.
            % save(fullfile(metrics_dir, name), 'cmos_all_data', '-v7.3');
            continue;
        end

        % ===== STAGE: feature / arrhythmia extraction =========================
        % Replace this block with your real extraction call(s), e.g. the
        % logic inside extract_arr_ps_dynamics.  After it runs, cmos_all_data
        % should carry df_map, pixel_size, ps_dynamics, etc.
        %
        %   cmos_all_data = run_feature_extraction(cmos_all_data);   % <-- yours
        %
        % Then validate, embed, aggregate, and persist:
        r = validate_metrics(cmos_all_data, name, species);
        cmos_all_data.qc = qc_attach(existing_qc(cmos_all_data), r);
        recs = [recs; r]; %#ok<AGROW>

        % save(fullfile(metrics_dir, name), 'cmos_all_data', '-v7.3');
    end

    % ===== Batch report =====================================================
    out_csv = fullfile(metrics_dir, ['batch_qc_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);
    [summary, detail] = qc_summarize(recs, out_csv);

    qc = struct('records', recs, 'summary', summary, 'detail', detail);
end


% ---- helper: existing embedded QC log, or [] if the field is absent ----
function qc = existing_qc(cmos_all_data)
    if isfield(cmos_all_data, 'qc')
        qc = cmos_all_data.qc;
    else
        qc = [];
    end
end


% ---- helper: did a given recording FAIL at a given stage? ----
function tf = any_fail(recs, name, stage)
    tf = false;
    for k = 1:numel(recs)
        if recs(k).recording == string(name) && recs(k).stage == string(stage) ...
                && recs(k).status == "FAIL"
            tf = true; return;
        end
    end
end
