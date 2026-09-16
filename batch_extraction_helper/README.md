# Batch feature extraction (headless)

Re-runs the Feature Extraction pipeline on every conditioned recording under
an input parent folder (any depth of sub-folders) without the app UI, writes
the `*-metrics.mat` files the other CADENCE modules consume into the same
relative folders under an output parent folder, and builds one Excel workbook
of per-recording metric medians grouped by experiment and by ZT.

```
input_root/ZT2_8am/R2/05-Rat-R2-20250511-8AM-Baseline-200ms-3v.mat
        -> output_root/ZT2/R2/05-Rat-R2-20250511-8AM-Baseline-200ms-3v-metrics.mat
           (a folder named ZT<n>_anything is written as ZT<n>; pass 'NormalizeZT', false to mirror names exactly)
output_root/cadence_recording_medians.csv      one row per recording (appended after every file)
output_root/CADENCE_median_summary.xlsx        Recordings | By experiment | By ZT | one sheet per metric | README
output_root/CADENCE_median_summary_*.csv       CSV copies of the three tables
output_root/cadence_batch_log.txt              timestamped log
```

## Run

Generic form, any study: supply the input parent folder and the output parent
folder.

```matlab
addpath('batch_extraction_helper');
T = run_cadence_batch('/path/to/conditioned_root', '/path/to/metrics_root');
T = run_cadence_batch(in, out, 'Folders', {'ZT2_8am/R2'});   % one sub-folder (and everything below it)
T = run_cadence_batch(in, out, 'DryRun', true);              % list what would run
```

From a shell (no desktop needed):

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -nodisplay -batch "addpath('batch_extraction_helper'); run_cadence_batch('/in', '/out')"
```

`run_g1_rat_batch.m` is the same call pre-filled for the G1 rat study; edit
`TEST_FOLDERS` there (`{}` = whole study). The run is
resumable: recordings already listed in the medians CSV whose metrics file
exists are skipped, so a long run can be interrupted and restarted. Use
`'Resume', false` to force re-extraction. `'DryRun', true` lists what would be
processed. Files are written to a `.part` name and renamed on completion, so
an interrupted save never leaves a truncated `-metrics.mat`.

## From raw camera files: convert, condition, extract

`run_cadence_pipeline` runs all three modules headlessly, starting from the
raw camera files. Every folder under the raw parent folder (any depth) that
holds `.tif` volumes or SciMedia `.gsh/.gsd` pairs is one recording (the Data
Conversion app's rules: the files are its cameras in directory order, or, for
legacy folders where every file ends in a camera letter A-D, one recording per
letter-stripped name). Each recording is converted, conditioned with the
Signal Conditioning app's default settings and extracted in turn.

```
raw_root/ZT2_8am/R2/05-Rat-R2-20250511-8AM-Baseline-200ms-3v/01_voltage.gsh  (+.gsd, 02_calcium.gsh ...)
        -> processed_root/converted/ZT2_8am/R2/05-Rat-R2-20250511-8AM-Baseline-200ms-3v.mat
        -> processed_root/conditioned/ZT2_8am/R2/05-Rat-R2-20250511-8AM-Baseline-200ms-3v-conditioned.mat
        -> metrics_root/ZT2/R2/05-Rat-R2-20250511-8AM-Baseline-200ms-3v-metrics.mat
processed_root/cadence_conditioning_log.csv    one row per conditioned recording: format, cameras,
                                               frames, frame rate, pacing found, polarity decision per
                                               camera, conditioning QC flags, stage errors, time
metrics_root/...                               medians CSV, workbook, log, excluded list as above
```

```matlab
addpath('batch_extraction_helper');
T = run_cadence_pipeline('/path/to/raw_root', '/path/to/processed_root', '/path/to/metrics_root');
T = run_cadence_pipeline(raw, proc, met, 'Folders', {'ZT2_8am/R2'});   % one sub-folder first
T = run_cadence_pipeline(raw, proc, met, 'DryRun', true);              % list recordings and what exists
```

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -nodisplay -batch "addpath('batch_extraction_helper'); run_cadence_pipeline('/raw', '/processed', '/metrics')"
```

`cadence_pipeline_ui` opens a window for the same run, laid out like the
other CADENCE modules: SELECT buttons for the raw, processed and metrics
folders, a list of the raw folder's sub-folders to restrict the run (none
selected = whole tree), the Signal Conditioning dropdowns with their defaults,
the field of view, the run options (dry run, resume, re-condition,
re-convert, keep converted files, condition arrhythmia-tagged recordings), a
RUN button, a STOP button that ends the run after the current recording, and
a console that echoes every batch log line. If the processed or metrics
folder is left empty, `<raw folder>_processed` and `<raw folder>_metrics`
next to the raw folder are used.

Conditioning reproduces the app's EXECUTE SIGNAL CONDITIONING with its
dropdown defaults, in the app's stage order: SNR map, drift correction, SVD
denoising (rank 8, 5 x 5 binning fallback), temporal filter [0, 50] Hz,
normalization, polarity check (post-stimulus deflection, shape fallback),
ensemble averaging on `analog1` (auto-detected beats when there is no pacing
channel), SNR >= 2 tissue mask, `conditioned` stamp and QC log. Spatial
binning and motion correction are off, as in the app. Change any of these with
`'ConditionOpts'` (fields `Drift, SVD, SVDRank, Binning, BinSize, FilterHz,
Motion, Normalize, Ensemble, MaskFloor`, see `cadence_condition_data.m`).
Extraction uses the per-camera adaptive SNR mask and the ensemble beat, as
described below; `'ExtractOpts'` passes through (FOV_mm etc.). `.tif`
recordings carry no frame rate; `'SamplingHz'` supplies it (default 1000 Hz,
the app's field default). `.gsh` headers always carry it.

The run is resumable at every stage: a recording already in the medians CSV
with its metrics file is skipped; one with a conditioned file is only
re-extracted; one with a converted file is only re-conditioned. `'Resume',
false` re-extracts everything from the conditioned files, `'Recondition',
true` rebuilds the conditioned files from the converted ones (use this after
changing `ConditionOpts`; a conditioned file is never conditioned again),
`'Reconvert', true` starts over from the raw files. A partial re-run
(`'Folders'` with any of these) replaces its own rows in the CSVs and keeps
the others. `'SaveConverted', false` skips writing the converted files if
disk is short; re-conditioning then reads the raw files again.

Arrhythmia-tagged recordings (see the next section) are converted and
conditioned, so they are ready for the arrhythmia-dynamics stages, but not
extracted; `'ConditionExcluded', false` skips them entirely. A conditioning
stage that throws fails the recording (nothing saved, `FATAL` in its row;
`'StrictConditioning', false` saves and extracts anyway, as the app would).

The conditioned tree is exactly what `run_cadence_batch(<processed_root>/conditioned,
<metrics_root>)` re-extracts, so the last stage can be re-run on its own.

## Cardiac maps as PDF

`'SaveMapsPDF', true` (a checkbox in the window) also writes
`<metrics_root>/<rel>/<name>-maps.pdf` for every extracted recording: the map
windows the Feature Extraction app pops up and saves after EXTRACT FEATURES,
one page per feature and camera in the app's order (activation, APD,
repolarization, AP rise, CV with propagation arrows, V-Ca delay, CaTD, Ca
decay, Ca rise, Ca tau, the batch-only APD50 / CaTD50, DF / RI / OI,
alternans substrate panels), each with the map over the camera image on the
app's robust colour scale, histogram, representative trace and statistics.
Pages are rendered headlessly (raster, like the app's PDF); the option works
with `run_cadence_pipeline`, `run_cadence_batch` and
`cadence_recompute_medians`, and the PDF path is recorded in the `Maps_pdf`
column of the medians CSV. On a resumed run, recordings already extracted but
without a PDF get one rendered from their saved metrics file, without
re-extraction; `cadence_recompute_medians(metrics_root, 'SaveMapsPDF', true)`
does the same for a whole existing metrics tree (about 15 s per file to load
plus a few seconds to render). `cadence_maps_pdf(metrics_struct, pdf_file)`
renders one recording.

## Excluded recordings

Recordings whose file name contains an arrhythmia tag are not extracted
(the pattern tolerates the misspellings present in the corpus, such as
"Arhhythmia", while not matching "alternans" or "EAD"). They call for the
arrhythmia-dynamics stages, not the paced metrics, and are listed in
`cadence_excluded_recordings.csv` with ZT, experiment, condition and CL, and the
onset sheet uses that list for the arrhythmia CL. Change or disable with
`'ExcludePattern'` (both `cadence_batch_extract` and `cadence_recompute_medians`).

## Analysis sheets

Besides the medians, the workbook carries three derived sheets (also written
as CSV):

- **Alternans onset** (one row per ZT × experiment × condition): the longest
  CL at which voltage and calcium alternans appear (`Onset_CL_V`,
  `Onset_CL_Ca`), which came first, the concordant → discordant → arrhythmia
  transition with the CL of each step, the functional refractory period
  (shortest CL still at 1:1 capture), and the same onsets read from the
  file-name tags (`Alternans`, `Dis`, `Arrhythmia`) for comparison. Criteria:
  alternans present when ≥ 10 % of tissue pixels have ratio > 0.10 and
  p < 0.05 and this also holds at the next shorter CL; discordant when the
  concordance ratio < 0.80; arrhythmia from the file tag, or from the metrics
  when the voltage DF exceeds the pacing rate by more than 15 % (the tissue is
  no longer driven by the stimulus; ratios below 1 are block, not arrhythmia),
  or when capture is lost with OI < 0.5. The per-recording sheet names the
  criterion that fired. Thresholds are options of `cadence_alternans_onset`.
- **Alternans per recording**: the classification behind it, one row per
  recording with its state string.
- **Circadian cosinor**: for every metric (per-experiment value at a reference
  CL, default 150 ms Baseline) and every onset measure, a 24 h cosinor fit
  (mesor, amplitude, acrophase, R², F-test p). Change the reference CL with
  `cadence_build_summary(..., 'RefCL', 125)`.

`cadence_recompute_medians(output_root)` rebuilds the CSV and workbook from
saved metrics files without re-extracting (about 15 s per file), for use after
editing the metric list or the analysis rules.

## What is computed, and how it matches the app

`cadence_extract_features` is a line-for-line port of the app's private
extractors (`extract_apd`, `extract_rise_time`, `extract_tau`,
`extract_v_c_delay`, `extract_cv`, `extract_act_time`, `extract_alternans`,
`extract_arr_complexity`, `extract_combo_masks`, `despeckle_act_map`) driving
the same shared helpers (`compute_lat_50`, `conduction_velocity`,
`analyzeAPAlternans`, `analyzeCaTransientAlternans`,
`assess_arrhythmia_substrate`, `cardiacSpectralMetrics`, `create_snr_mask`,
`extract_beat_frames`, `normalize_data`). Output cell layouts
(`ep_metrics.<field>{camera, slot}`) are identical to the app's, and
`window`, `data_masks`, `FOV` are written as `compile_data` writes them.

Settings reproduce the app's batch defaults on the lab's dual-camera rig
(these are also what the existing metrics files on the drive were made with):

| Setting | Value |
|---|---|
| Analysis window | ensemble-averaged beat (`CAM<n>_average`); no window drawn |
| Mask | adaptive SNR mask per camera, `create_snr_mask(snr, [], true)` |
| CAM1 (voltage) | LAT (50 % upstroke), APD80, repolarization 80 %, rise 20–90 %, local CV (inverse gradient, FOV 20 mm) |
| CAM2 (calcium) | V–Ca delay, CaTD80, decay time 80 %, rise 20–90 %, tau fit 90 %→20 % plus model-free 1/e |
| Alternans | AP alternans on CAM1, Ca alternans on CAM2, combined substrate; beats from `analog1` |
| DF / RI / OI | CAM1 (voltage) over the full record |
| Representative pixel | image centre |

Change any of these through `ExtractOpts` (see the header of
`cadence_extract_features.m`). One deliberate addition beyond the app's
output: APD50 and CaTD50 maps are stored as `ep_metrics.apd50_data` and
`ep_metrics.ca50_data` (same 6-slot layout, same engine as APD80; option
`ExtraAPDLevels`, pass `[]` to omit). The app's own fields are unchanged. Arrhythmia-dynamics stages (wavefront and
phase-singularity tracking) are not run; they are not part of the summary.

Parity with the app was checked two ways: (1) `mv_reextract_check` between an
app-produced metrics file and the headless output of the same recording gives
bit-identical activation, APD, repolarization, CaTD and Ca decay-time maps
and identical CV medians (the rise-time, tau, V–Ca delay and alternans maps
differ only because those algorithms were rewritten after the app file was
made); (2) a line-level diff of each ported extractor against the app source
shows the computational lines are identical. To re-check on a recording
extracted by the current app build, run BATCH PROCESSING in the app on that
one file into a scratch folder, then:

```matlab
mv_reextract_check('/scratch/<name>-metrics.mat', '/metrics_root/ZT2/R2/<name>-metrics.mat')
```

Every field should report "identical". If the app's extractors are edited,
port the change here too and repeat.

## Spreadsheet columns

Every value is the median over finite pixels of the masked map for that
recording. Definitions are in `cadence_metric_labels.m` and on the README
sheet of the workbook. Conduction velocity is pooled over the tissue interior
(8 px erosion), the same gate Signal Analysis and the Conduction Velocity
module apply; two CV columns are given. Both are the inverse-gradient method
(v = ∇T/‖∇T‖²) and differ only in how ∇T is estimated: the map saved by Feature
Extraction uses finite differences, the recomputation uses a local polynomial
gradient fit (polyfit), the CV module's default estimator.

Column labels reused from the previous Python summary (`extract_medians.py`)
are unchanged: Activation time (V) [ms], APD (V) [ms], Repolarization time (V)
[ms], AP rise time (V) [ms], Voltage-Ca delay [ms], Ca transient duration
[ms], Ca decay time [ms], Ca rise time [ms], Ca decay tau [ms], Dominant
frequency (V) [Hz], Regularity index (V), Organization index (V), APD spatial
gradient (V), APD restitution slope (V), APD alternans ratio (V), AP amplitude
ratio (V), Alternans phase (V).

Renamed or replaced: the old "(Ca)" copies of the substrate columns (APD
spatial gradient, restitution slope, APD alternans ratio, Alternans risk) were
duplicates of the voltage values and are replaced by true calcium alternans
metrics (CaT amplitude alternans ratio, CaTD50/80 alternans ratio, Ca release
and load alternans, Ca spectral index, Ca–AP coupling). "Alternans risk (V/Ca)"
is now a single "Arrhythmia risk (composite)" from the combined Vm+Ca
substrate. "Conduction dv/dt ratio" is now "dV/dt alternans ratio (V)".

Metadata columns come from the folder and file names: the nearest ancestor
folder named `ZT<n>...` → ZT; the nearest ancestor folder named `R<n>` →
Experiment (else the immediate parent folder);
`<run>-Rat-R<rat>-<date>-<clock>-<Baseline|IP>-<CL>ms-<stim>v[-<tag>]`.

## Beat-windowed columns and ensemble validity

The timing metrics above are measured on the ensemble-averaged beat. At fast
pacing the calcium transient may not relax within the cycle; the ensemble
window then references the "baseline" to the previous transient's tail and
CaTD, Ca rise, tau and V–Ca delay collapse (a negative V–Ca delay is the
tell). Rather than substituting a different estimator at those CLs, which
would confound method with pacing rate, every recording also gets a
**beat-windowed** column set taken from the per-beat maps the alternans
analysis already stores: each beat in its own stimulus-to-stimulus window,
referenced to its own trough. Ca time-to-peak, decay-50 and decay-80 from the
peak, their onset-referenced sums (CaTD50/80), the diastolic level, and
beat-windowed APD50/APD80 for voltage. Values are the per-pixel median over
beats, then the median over pixels. Use these consistently across all CLs.

At fast pacing the transient legitimately becomes more symmetric (rise and
decay speeds converge), so a polarity detector run on the ensemble beat drifts
toward zero and starts returning arbitrary signs. Do not use one to
"auto-correct" fast-paced recordings; the sign carries no information there.

Two flags, `Ca ensemble metrics valid` and `V ensemble metrics valid`, mark
recordings where the ensemble values pass every test: the same guard the app
applies to the averaged beat (`check_analysis_window`, in
`feature_extraction_helper` — window not wrapped, decay fraction at least 0.5),
plus outcome checks the guard cannot make because it never sees the extracted
values (non-negative V–Ca delay, enough valid pixels, plausible duration,
within 50 % of the beat-windowed value). Rows flagged 0 should have their
ensemble Ca (or V) timing cells treated as missing and the beat-windowed
columns used instead.

Note the difference in behaviour between the two tools, which is deliberate.
The app **refuses** to write a feature whose window fails the guard, because
its user is looking at one recording and a missing map is an unmissable
signal. This batch tool **extracts anyway and flags**, because a study-wide
spreadsheet is more useful when the questionable value is visible next to the
flag and the beat-windowed alternative. Re-run `cadence_recompute_medians` to
refresh the flags after changing the guard; it does not re-extract.

## QC columns

Three columns exist only to catch the failure modes that matter for the
downstream questions:

- **V and Ca maps identical (QC flag)** — 1 if the APD80 map equals the CaTD80
  map, the two activation maps are equal, or the two camera stacks are the
  same array. Any 1 means the calcium results were not read from the calcium
  camera and the row should not be used.
- **Pacing rate measured [Hz]** — from the stimulus channel (`analog1`).
- **Capture ratio (DF / pacing rate)** — voltage dominant frequency over the
  pacing rate: ≈1 is 1:1 capture, ≈0.5 is 2:1 block. Within one experiment and
  condition, the shortest CL whose capture ratio is still ≈1 is the functional
  refractory period; the first CL below it that drops to ≈0.5 marks loss of
  capture.

Note that the restitution slope, APD gradient and concordance columns are
voltage (APD vs diastolic interval) quantities by definition; CADENCE computes
no calcium restitution. The calcium alternans columns (CaT amplitude, CaTD50/80
ratio, release/load) come from `analyzeCaTransientAlternans` on CAM2 and are
independent of the voltage ones; Ca–AP coupling r and the in-phase fraction are
the cross-camera measures.
