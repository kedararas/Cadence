# Validation study tooling

Two tracks, with different ground truths:

- **`mv_mark` and friends** — blinded, reviewer-driven manual marking, for metrics
  where a human eye is the only available reference (APD, CaTD, alternans). This is
  the approach KairoSight-3.0 used, with the blinding and independence made explicit
  rather than left unstated.
- **`mv_paced_df`** — fully automated validation of dominant frequency against the
  pacing rate recorded in the stimulus channel. No human in the loop, and it scales
  to every paced recording you own.

The folder name says "manual" for historical reasons; it holds both.

## `mv_paced_df` — automated DF validation

Under S1S1 pacing, DF must equal the pacing frequency. The stimulus channel carries
that truth in every file, so this is exact, external, and free.

```matlab
R = mv_paced_df('/Users/kedararas/Desktop/Test Data', 'Stride', 8, 'Output', 'df.csv');
```

Files are 2–8 GB, so the camera stack is read as a strided HDF5 slab (~33 MB at
stride 8) and the spatial-mean trace is cached in a small sidecar. The first pass
costs a few seconds per file; every re-run afterwards is instant.

Run it **with the shipped default band first** — that is what a user actually gets,
and it is the primary result. A fixed wide band (`'Band', [0.5 50]`) is a legitimate
secondary run because it is identical for every recording. Never set a band around a
specific recording's known pacing rate; that is circular.

### What the first pilot found

On 33 recordings (10 rat, 23 human):

- With a band that reaches the rhythm, the DF engine picked the **correct FFT bin in
  every single recording**, across 1–14.3 Hz. All returned values are exact multiples
  of the 0.244 Hz bin, and the maximum deviation from the pacing rate is one half-bin
  — the residual is frequency quantization, not algorithmic error.
- **`cardiacSpectralMetrics` cannot measure any rhythm below 2 Hz.**
  `scanLo = 3` and `lo = max(2, 0.5*f0)` (lines 176–188) put a hard floor on the
  search band, so human and large-animal recordings paced at 1000 ms / 60 bpm return
  a plausible-looking wrong answer (2.2–2.9 Hz instead of 1.0 Hz) rather than failing
  loudly. RI collapses to 0.05–0.2 when this happens, so RI is a usable detector.
  The adaptive band is intentional and correct for VF; the fix is to lower the floor,
  not to revert to fixed bands.
- Reporting agreement as "% within 1%" is meaningless below ~5 Hz, because 1% of 1 Hz
  is finer than the 0.244 Hz bin. Report absolute error against the bin width instead.
  Sub-bin parabolic interpolation of the spectral peak would cut the error by roughly
  an order of magnitude and is a few lines.

---

## Manual marking harness

## Why this is a separate folder

The harness **shares no code with `feature_extraction_helper`**. It does not import
the activation-time detector, the repolarization threshold logic, or anything else
under test. It loads `.mat` files with plain `load()` and does its own field checks
rather than calling `load_cmos`, and it computes its own percentiles rather than
using the Statistics Toolbox.

That duplication is deliberate. A reference that calls the code it is meant to check
is not a reference, and this is the first thing a reviewer will look for. Keeping it
visibly standalone costs about a hundred lines and settles the question at a glance.

## Workflow

```matlab
addpath('manual_validation_helper');

% 1. Draw the CANDIDATES — ONCE per recording. Twice the pool size, from
%    pixels above the same adaptive SNR floor Feature Extraction masks with.
%    Reviewers mark CAM1_average by default, the same ensemble-averaged beat
%    feature extraction uses. 'grid' lays the pixels on a lattice so they can
%    be interpolated into a manual map; 'stratified' draws them by SNR tercile.
manifest = mv_sample_pixels('rat-conditioned.mat', 'CAM1', 100, 'Sampling', 'grid');

% 2. The LEAD reviewer builds the pool: candidates are shown in manifest
%    order, unmarkable ones are skipped, and the session stops by itself at
%    100 accepted. Those 100 are the pool; the skips are recorded.
mv_mark(manifest.manifest_file, 'AA', 'Lead', true);

% 3. Every other reviewer marks the pool, independently and blinded
%    (refused until the lead session has finished)
mv_mark(manifest.manifest_file, 'SP');
mv_mark(manifest.manifest_file, 'EM');

% 4. Derive metric values from the saved clicks
tbl = mv_derive({'..._AA_marks.mat', '..._SP_marks.mat', '..._EM_marks.mat'});

% 5. Compare against CADENCE, with inter-observer spread as the benchmark
stats = mv_compare(tbl, 'rat-metrics.mat', 'Field', 'apd_data');
```

In the breadth arm the single reviewer *is* the lead, so steps 2 and 3 collapse into
one session. The marking window shows a progress line and bar throughout — for the
lead, accepted against the target plus the candidate count; for everyone else,
marked, skipped and remaining.

For calcium, point step 1 at the calcium camera and mark that channel separately:
`mv_sample_pixels('pig-conditioned.mat', 'CAM2', 100, 'Sampling', 'grid')`. Rabbit
recordings are voltage-only, so they have no calcium arm.

Run `mv_selftest` first. It builds a synthetic recording whose true APD80 is known
analytically, simulates three reviewers, runs the whole non-interactive path, and
asserts that a deliberately injected bias is recovered. If it passes, everything
downstream of the human clicking is known good.

### What runs per recording, and what runs once

Only the first group involves reviewers. The rest are automated or synthetic, run
**once for the whole study**, and are not repeated per recording — a common source of
confusion, since they live in the same folder.

| | When | What it covers |
|---|---|---|
| `mv_sample_pixels` → `mv_mark` (lead) → `mv_mark` ×2 → `mv_derive` → `mv_compare` | once per **recording × channel** | activation, APD/CaTD at every level, amplitude |
| the same chain with `'Mode','alternans'` | once per **alternans-positive recording × channel** | amplitude alternans |
| `mv_paced_df` | **once**, across every paced recording you own | dominant frequency, against the pacing hardware |
| `mv_synth_cv`, `mv_cv_envelope`, `mv_lat_quantization` | **once** | conduction velocity, against a synthetic planar wave |
| `mv_rise_test` | **once** (and after any edit to `extract_rise_time`) | rise-time algorithm, against analytic truth |
| `mv_lat_test` | **once** (and after any edit to `compute_lat_50`) | activation-time algorithm, against analytic truth |
| `mv_selftest` | **once**, before recruiting reviewers | the harness itself |

`mv_rise_test` and `mv_selftest` are regression tests, not study instruments: they
take no reviewer input and produce no per-recording result. `mv_rise_test` earns its
place in the coverage table the same way `mv_synth_cv` does — as synthetic ground
truth for a metric no human can reference — but it is a single pass/fail run, not
something to repeat per heart.

## Mark what the software measured: the ensemble-averaged beat

`mv_sample_pixels` defaults to `'Source', 'ensemble'`, which points the reviewers at
`CAM<n>_average` rather than the raw stack.

This follows the pipeline. Ensemble averaging is CADENCE's recommended conditioning
step, and feature extraction takes `avg_cmos_data` whenever that checkbox is ticked —
so every metric except alternans is derived from the averaged beat. Marking the same
array is what makes the comparison like-for-like. It also removes beat-to-beat
variation as a source of apparent disagreement, which would otherwise inflate the
limits of agreement without telling you anything about the software.

**And it dissolves a problem that would otherwise block half the metrics.** CADENCE
stores some metrics as durations and others as times measured from the start of its own
analysis window:

- **Durations** — `apd_data`, `ca_data`, `ap_rise_times`, `ca_rise_times`, `ca_tau`,
  `vc_delay`. Differences of two times, so any origin cancels on both sides.
- **Window-relative times** — `act_times`, `rep_data`, `ca_rep_data`.

On the raw arm, a reviewer clicks in absolute recording time while the software counts
from its own window start, so differencing them uncorrected gives a constant bias of
`(2 - start_frame)` frames — for a window starting at frame 1001 at 1 kHz, −999 ms on
every pixel, with immaculate-looking scatter and correlation. A silent failure.

On the ensemble arm that offset does not exist: both sides start at frame 1 of the same
averaged array, and the only residual is the 1-based index convention (`compute_lat_50`
returns a 1-based fractional frame; `mv_mark`'s axis is `(frame-1)*dt`), which is a
single frame and is corrected automatically. `mv_selftest` verifies activation time
validates to a bias of +0.05 ms with no window information involved at all.

`mv_compare` classifies every field, applies the right correction for the arm in use,
and **refuses** rather than guessing when it cannot verify one. It also sanity-checks
the corrected values against the marked window and warns if most fall outside — the
metrics file does not record which array feature extraction actually used
(`ensemble || isempty(oap_window)` decides at run time and nothing is saved), so that
check is the only available guard against the two sides having looked at different data.

**Alternans is the exception, and is refused on this arm.** Averaging is precisely what
removes the beat-to-beat alternation, so `'Source', 'ensemble'` with `NumBeats > 1`
errors. Run the alternans session with `'Source', 'raw'`.

```matlab
% the workhorse arm — ensemble by default
manifest = mv_sample_pixels('rat-conditioned.mat', 'CAM1', 100, 'Sampling', 'grid');

% the alternans arm — raw, two beats
alt = mv_sample_pixels('rat-conditioned.mat', 'CAM1', 75, 'Source', 'raw', ...
                       'NumBeats', 2, 'MetricsFile', 'rat-metrics.mat');
```

**One thing to state carefully in the methods.** `CAM<n>_SNR` is computed on the RAW
stack, before averaging (Signal Conditioning runs it first). So under this arm the SNR
strata describe single-beat signal quality while the reviewers are marking a much
cleaner averaged trace. That is the useful reading — agreement as a function of the
underlying data quality — but it must be reported that way, not as the SNR of the trace
on screen.

For the raw arm, pass `'MetricsFile'` so the reviewers mark the same beat CADENCE
analysed; it reads only the window, never a metric value, so blinding is unaffected.

## Addressing the right map: `Camera` and `Slot`

`ep_metrics` entries are `{camera, slot}` — the row is the camera (the `extract_*`
drivers loop `for i = 1:num_files` and pass `i` straight in as the camera argument),
and the column is a data slot: 1 is the masked map, 2 the background image, 5 the
camera label, 6 the unmasked map.

The old `FileIndex`/`CamIndex` names described this backwards and are now rejected with
an error. They were dangerous rather than merely confusing: `CamIndex, 2` returned the
**background image**, which is numeric and correctly sized, so it sailed through and
produced confident nonsense. The old auto-resolve-by-size had the same failure — it
took the first `[nr nc]` match, which on a rig with two voltage cameras is camera 1's
map regardless of who was marked.

`mv_compare` now takes the camera from the manifest (the camera the reviewers actually
marked) and addresses the cell directly, erroring if the slot is empty or the wrong
shape. Override with `'Camera'` and `'Slot'`.

Two smaller corrections in the same pass. `act_times` is masked with **zeros** rather
than NaN, so off-tissue pixels arrived as a finite 0, passed the `isfinite` filter and
dragged the bias negative; those are now treated as missing. And the inter-observer
benchmark is computed on **all** marked rows rather than the software-resolved subset —
conditioning the human noise floor on the pixels CADENCE handled removes exactly the
hard ones from the human side and flatters the comparison the report leads with.
`mv_derive` also now refuses to pool marks files from different manifests, since the
`pixel` column is an index into the manifest and pooling makes it mean different
locations for different reviewers.

## The four design decisions that matter

**Blinding.** `mv_mark` never opens a `*-metrics.mat` and never displays a CADENCE
value — not even as a marker "for reference". Reviewers anchor to a visible answer,
and agreement measured that way is manufactured. Row and column are withheld too, so
a reviewer who has already seen the APD map cannot place the pixel on it from memory.
Comparison happens in `mv_compare`, afterwards.

**Human supplies judgement, code supplies arithmetic.** The reviewer's click fixes a
*time*; the harness reads the trace value there. It never searches for a nearby local
maximum, because that search *is* the automated detection under test. The single
exception is the diastolic baseline, taken as a median over ±5 ms, because one sample
is noise-dominated — and that is recorded in the output so it can be reported.

For APD the reviewer marks baseline and peak, the harness draws a guide line at the
repolarization level, and the reviewer then marks where the trace crosses it. This
converts an unreliable judgement ("where is 80% repolarized?") into a reliable one
("where does the trace cross this line?") without automating the decision.

**Clicks out of order are slips, and are handled by the clicks alone.** The 50% guide
line is crossed twice, so an activation click can land on the downstroke. That gives a
short but positive APD, which the APD <= 0 check does not catch. In the first study
this happened 4 times in 3,600 marks, all by one reviewer. `mv_mark` now refuses an
activation click after the reviewer's own peak click, or a repolarization click before
it. `mv_derive` drops such marks from older files and lists them in
`tbl.Properties.UserData.dropped`. Neither rule looks at a CADENCE value, so applying them
cannot favour the software. The rule was written after those slips were seen, so report
results with and without it.

**Activation gets a guide line too, and it must.** CADENCE times activation at the
50% baseline-to-peak crossing of the upstroke (`compute_lat_50`), not at maximum
dV/dt. Asking a reviewer for "the upstroke" invites them to mark max dV/dt, and the
offset between the two definitions would then be booked as software-versus-human
disagreement when it is really the two sides measuring different quantities. So the
third click is marked against a 50% guide line, and `mv_mark`'s `ActPercent` defaults
to 50 to match the code under test. This matters twice over: APD is measured from the
activation click, so it inherits any error there.

**Candidates, then a pool — the lead decides markability, not the code.** A pixel
can clear an SNR floor and still be unmarkable: a distorted morphology, a motion
artefact, a double hump. Left as is, reviewers skip such pixels and the study ends up
short of its count, unevenly across reviewers. The harness must not pre-screen them
itself — "is this a usable action potential" is exactly the judgement under test —
so `mv_sample_pixels` draws **twice** the pool size and leaves the pool *open*. The
**lead reviewer** then runs `mv_mark(..., 'Lead', true)`: candidates come up in
manifest order, the lead marks or skips each, and the session stops on its own when
the target is reached. Accepted pixels become the pool (`manifest.active`), the
lead's skips are recorded and reportable as a skip rate, and candidates never reached
are discarded. Other reviewers, `mv_derive` and `mv_compare` see the pool only, and a
non-lead session is refused while the pool is open. The lead is blinded exactly like
everyone else. Candidate order is built so the pool is sound wherever the lead
stops: `'stratified'` interleaves strata round-robin and the lead fills a
**per-stratum quota** (a dim stratum is skipped more, and a global count would let
the bright ones crowd it out); `'grid'` presents the base lattice first, then a
lattice offset by half a stride (the cell centres), so the pool stays evenly spread.
Grid strata are assigned by rank once the pool is final.

The population this validates is therefore *pixels a trained observer can mark
above the analysis-mask floor*, and the Methods must say so. The skip rate is part
of the result.

**Candidates are drawn above the adaptive SNR floor.** Feature Extraction masks its
maps with `create_snr_mask(snr, [], true)` — floor = max(1.5, half the median tissue
SNR) — so pixels below it have no CADENCE value and were being dropped in
`mv_compare` anyway. `mv_sample_pixels` now applies the same rule (re-implemented in
one line, so the harness still imports nothing from the pipeline) and records the
resolved floor in the manifest. Consequence for the SNR envelope: the lowest stratum
means "usable but dim", not "noise". `'SNRFloor', 0` restores the old SNR > 0 rule.

**Seeded sampling, in one of two layouts.** Either way the seed makes the pixel list
a pre-registered choice — hand-picking pixels, or redrawing until the numbers look
good, invalidates the exercise.

- `'stratified'` (default) draws equal numbers from each SNR tercile. This lets
  `mv_compare` report agreement *as a function of SNR*, an operating envelope users
  actually need, and pre-empts the obvious criticism that validation used only clean
  signals.
- `'grid'` places the same number of pixels on a regular lattice instead. It costs
  the same marking effort but the marks can be **interpolated into a manual map**,
  which scattered pixels cannot be — so the same clicks yield the map panels, the
  agreement statistics *and* the SNR envelope. SNR strata are assigned afterwards
  (post-stratified, by rank so ties cannot empty a stratum). Points that fall off
  the tissue mask are dropped rather than nudged to a neighbour, since moving them
  would bend the lattice and cost the even spacing that is the whole point.

  The honesty rule that goes with grid mode: interpolate for the **map panels only**.
  Every statistic is computed at the marked points, and the difference map is plotted
  as dots at those points, never interpolated.

  A square lattice can only yield certain counts (8×8, 9×9, 11×11), so `n_pixels` is
  a target rather than a quota. The stride chosen and the count actually drawn are
  printed and stored in the manifest.

**Raw clicks are the record.** `mv_mark` saves click coordinates, not metrics.
`mv_derive` computes values from them afterwards, so a definition can change —
APD80 to APD50, a different alternans normalisation — without anyone re-marking.

## Reading the output

`mv_compare` reports software-versus-manual agreement *and* inter-observer spread.
The second is what makes the first interpretable: if software-versus-human
disagreement is no larger than human-versus-human disagreement, the software sits
inside the noise floor of expert manual measurement. That claim does not require
anyone to accept that 5 ms is the right tolerance, which is why it is the one to lead
with.

## Scope

Manual marking is the right ground truth for repolarization time, APD and CaTD at any
percentage, amplitude, and amplitude-alternans ratios. It is **not** the right tool
for everything. This is the coverage table:

| Metric | Ground truth | Notes |
|---|---|---|
| Activation time | This harness | Marked against a 50% upstroke guide, matching `compute_lat_50`. Not weak, provided the guide is used — see above. |
| APD 30/50/80/90 | This harness | One marking pass covers every level: `mv_derive` recomputes from raw clicks, so no one re-marks for a new definition. Repolarization time and APD are **not independent evidence** — same four clicks. |
| Repolarization time | This harness | Comes from `extract_apd` (which writes `rep_data` from its own columns 7/8/9), so it is min-max normalized and sub-frame interpolated — clean. The legacy `extract_rep` is a fallback only, reachable when APD is not selected; see below. |
| CaTD at any level | This harness | Voltage and calcium are separate marking passes on separate cameras. |
| Amplitude alternans | This harness | `'alternans'` mode. Validate on recordings that actually alternate — near zero the ratio is noise-dominated and Bland–Altman degenerates. |
| APD / CaTD alternans | Not yet | Needs a 2-beat APD mode in `mv_mark`. Note the software averages over *all* beat pairs while a reviewer sees two, so the estimators differ. |
| dV/dt alternans, spectral alternans index | None | SAI is an FFT over the beat series; there is nothing for a human to mark. |
| AP / Ca rise time | Synthetic (`mv_rise_test`) | **Not** manually marked — `mv_mark` has no rise-time mode. Covered by analytic ground truth instead, the same way CV is. Valid only from the 2026-09-01 fix below; earlier builds carry a ~+1 frame bias and fail at low SNR. Adding a manual arm would need a new mode, and is only worth it for calcium — see below. |
| Decay constant (τ) | None | An exponential fit, not a crossing. A reviewer cannot supply a reference for it. |
| Conduction velocity | Synthetic planar wave | `mv_synth_cv`, `mv_cv_envelope`, `mv_lat_quantization`. Validated — just not manually. |
| Dominant frequency | Paced recordings | `mv_paced_df`. DF must equal 1/CL, an exact truth already in your data, at n in the thousands. |
| Rotor / phase singularity | None available | Report as qualitative demonstration, or build an analytic rotating spiral with a known PS track. |

`ap_waveform` inside `mv_selftest` is the AP model to reuse for the synthetic-wave
work; propagate it across a grid at a known velocity and the CV ground truth follows.



## Study design: how to spend the marking hours

The binding constraint is reviewer time, and the unit of work is
**recording × channel × marking mode**, not recording. Running every mode on every
recording multiplies out to roughly 6,300 pixel visits per reviewer (~20 hours), which
nobody finishes — and reviewers who quit halfway leave you a biased subset. The
allocation below buys the same claims for about a quarter of that.

**Two arms, because they buy different things.**

- **Overlap arm — 4 recordings (one per species), all 3 reviewers.** Its only job is
  to estimate one variance component: how much two humans disagree on the same pixel.
  That saturates fast. What it needs is *species* coverage, not recording count, since
  5 ms is ~10% of a rat APD80 but ~2.4% of a pig's. A fifth or twelfth overlap
  recording just tightens a CI on a number used qualitatively.
- **Breadth arm — 12 recordings, one reviewer each.** Its job is to sample the
  replicate unit, and the replicate unit is the **heart**. Pixels within a heart are
  spatially correlated pseudo-replicates, so the marginal information in the 150th
  pixel of a heart already marked is close to nothing while the first pixel of a new
  heart is worth a lot. A breadth recording costs a third of an overlap recording.

Sixteen recordings total is also what the corpus allows: rabbit has exactly four.

**Rotate reviewers across species in the breadth arm.** Give each reviewer one
recording per species, not four of one species — otherwise reviewer identity is
confounded with species and any species difference in agreement is uninterpretable.

**~100 pixels per recording, on a grid.** Cheaper than 150 and strictly more useful,
for the spatial-correlation reason above plus the map panels. The lead visits more
than 100 — 100 plus whatever they skip — so the lead's budget grows by the skip rate;
everyone else marks exactly the pool.

**Alternans is a separate, targeted session.** It is a conditional metric: it only
exists where the tissue alternates. Run it on the two or three recordings that
actually do, all three reviewers, ~75 pixels — under an hour each, and it yields an
inter-observer benchmark for alternans as well.

| Arm | Recordings | Channel-sessions per reviewer | Pixels |
|---|---|---|---|
| Overlap (3 reviewers) | 4 | 7 (4 voltage + 3 calcium) | 700 |
| Breadth (1 reviewer each) | 12 | 7 (4 voltage + 3 calcium) | 700 |
| Alternans (3 reviewers) | 2–3 | ~3 | 225 |

About 1,625 pixel visits per reviewer, mostly four-click: 5–6 hours. Calcium counts
are 3-of-4 and 9-of-12 because rabbit is voltage-only.

The cost of the split: the breadth arm has no within-recording inter-observer
estimate, so a disagreeing breadth recording cannot be separated from its reviewer.
Rotation removes the systematic part; the rest is a stated limitation.


