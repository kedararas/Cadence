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

% 1. Draw the pixel sample — ONCE per recording, shared by all reviewers
manifest = mv_sample_pixels('rat-conditioned.mat', 'CAM1', 150, 'BeatIndex', 5);

% 2. Each reviewer marks the same pixels, independently and blinded
mv_mark(manifest.manifest_file, 'AA');
mv_mark(manifest.manifest_file, 'SP');
mv_mark(manifest.manifest_file, 'EM');

% 3. Derive metric values from the saved clicks
tbl = mv_derive({'..._AA_marks.mat', '..._SP_marks.mat', '..._EM_marks.mat'});

% 4. Compare against CADENCE, with inter-observer spread as the benchmark
stats = mv_compare(tbl, 'rat-metrics.mat', 'Field', 'apd_data');
```

Run `mv_selftest` first. It builds a synthetic recording whose true APD80 is known
analytically, simulates three reviewers, runs the whole non-interactive path, and
asserts that a deliberately injected bias is recovered. If it passes, everything
downstream of the human clicking is known good.

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

**Stratified, seeded sampling.** Pixels are drawn from SNR terciles with a fixed seed.
Stratification lets `mv_compare` report agreement *as a function of SNR*, which is an
operating envelope users actually need and which pre-empts the obvious criticism that
validation used only clean signals. The seed makes the pixel list a pre-registered
choice — hand-picking pixels, or redrawing until the numbers look good, invalidates
the exercise.

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
percentage, amplitude, and alternans ratios. It is **not** the right tool for
everything:

| Metric | Ground truth to use |
|---|---|
| APD, CaTD, repolarization, alternans | This harness |
| Activation time | This harness, but weak — humans judge max dV/dt poorly. Prefer the synthetic wave. |
| Conduction velocity | Synthetic planar wave at known velocity; cross-engine agreement (gradient / polyfit / circle) |
| Dominant frequency | Paced recordings — DF must equal 1/CL, an exact truth already in your data |
| Rotor / phase singularity | No ground truth available. Report as qualitative demonstration. |

`ap_waveform` inside `mv_selftest` is the AP model to reuse for the synthetic-wave
work; propagate it across a grid at a known velocity and the CV ground truth follows.

## Draft Methods text

> Automated feature extraction was validated against blinded manual measurement.
> For each recording, N pixels were drawn at random from within the tissue mask,
> stratified into equal groups by signal-to-noise ratio using a fixed random seed
> fixed before marking began. Three reviewers independently marked fiducial points on
> each sampled pixel's conditioned optical trace — diastolic baseline, peak,
> activation, and the crossing at X% repolarization — using a purpose-built marking
> tool that shares no code with the CADENCE feature-extraction routines and that
> displayed neither the automated result nor the pixel's location. Pixels were
> presented in an order randomized per reviewer. Metric values were computed from the
> recorded fiducial times after marking was complete. Agreement between automated and
> manual values is reported as bias and 95% limits of agreement (Bland–Altman), the
> proportion of pixels within X ms, and Pearson correlation, overall and by SNR
> stratum. Pairwise inter-observer differences on the same pixels are reported as the
> benchmark against which automated agreement should be judged.

Fill in N, X, and the tolerances from your run. State the conditioned-not-raw choice
explicitly: this validates feature extraction, not signal conditioning, and mixing
the two makes any disagreement unattributable.
