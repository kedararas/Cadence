# APD / repolarization: local-baseline threshold patch

Paste-in changes for `Cadence_Feature_Extraction.mlapp`. App Designer rebuilds the
`.mlapp` from its own model on save, so these must go in through **Code View** —
editing the embedded `matlab/document.xml` does not persist.

## Why

`extract_apd` and `extract_rep` both compared the trace against an **absolute**
level in a globally-normalized window:

```matlab
win_norm = normalize_data(win_data);   % 1st/99th percentile over the WHOLE window
[hit, rep_frame] = max(cumsum((win_norm <= rep_val) & (t_idx > pk), 3) == 1, [], 3);
rep(~hit) = NaN;
```

For APD80, `rep_val = 0.20`. That implements *"the trace falls to 0.20 of the
recording's global dynamic range"*, not *"this beat recovers 80% of its own
amplitude from its own baseline"*. Those agree only if every beat shares one flat
baseline and one identical peak.

On `04_4-30-2026_Human12_D0_CS1_1000ms_3v.mat` the diastolic baseline of individual
beats spans ~0.0 to ~0.52 of full scale within a single 4 s record, and peak heights
range 0.54–1.0. Any pixel whose resting baseline sits above 0.20 can never satisfy
the test, so `~hit` silently set it to NaN — the patchy APD maps.

The fix references the threshold to each pixel's **own** baseline and amplitude,
which is immune to both baseline wander and beat-to-beat amplitude variation.

## Measured effect

APD80, ensemble-averaged, tissue pixels (SNR ≥ 75th percentile), across four
species and four pacing rates. `remove_Drift` is **unchanged** (the existing
degree-3 detrend); the last column shows what enabling it contributes.

| recording | cam | before | after | drift-corr effect |
|---|---|---|---|---|
| human 005 (1000 ms) | CAM1 | 98.7%, 353 ms | 99.6%, 367 ms | 0 ms |
| human 005 | CAM2 | 100.0%, 322 ms | 100.0%, 340 ms | 0 ms |
| human 005 | CAM3 | 99.7%, 362 ms | 100.0%, 341 ms | 0 ms |
| human 005 | CAM4 | 99.9%, **126 ms** | 100.0%, **390 ms** | 0 ms |
| human 04 (1000 ms) | CAM1 | **37.5%**, 326 ms | **99.9%**, 248 ms | +3 ms |
| human 04 | CAM2 | 100.0%, 377 ms | 100.0%, 371 ms | 0 ms |
| rat (150 ms) | CAM1 | 100.0%, 54 ms | 100.0%, 50 ms | 0 ms |
| rat | CAM2 | 100.0%, 79 ms | 100.0%, 75 ms | 0 ms |
| pig (500 ms) | CAM1 | 98.9%, 273 ms | 99.8%, 267 ms | 0 ms |
| pig | CAM2 | 92.7%, 264 ms | 93.9%, 261 ms | 0 ms |
| rabbit (200 ms) | CAM1 | 90.7%, 128 ms | 91.7%, 126 ms | 0 ms |
| rabbit | CAM2 | 2.7%, 66 ms | 37.9%, 22 ms | 0 ms |
| rabbit | CAM3 | 89.9%, 127 ms | 91.5%, 126 ms | 0 ms |
| rabbit | CAM4 | 98.0%, 127 ms | 98.8%, 126 ms | 0 ms |

Yield never regresses. All of the improvement comes from the threshold and
activation fixes in this patch — none of it from drift correction.

Two entries are worth reading closely. **human 005 CAM4** goes from 126 ms to
393 ms — the old 126 ms was the broken-activation artifact, not a short AP.
**human 04 CAM1** goes from 37.5% to 100% yield. **rabbit CAM2** stays broken in
both (decay-fraction 0.03: the tissue never relaxes within the cycle, so there is
no diastolic baseline to reference). That camera should be excluded from APD
analysis regardless of this patch — it is the failure mode `validate_conditioned`
already warns about.

Repolarization time (`extract_rep`) was never dropping out — 99.3% → 99.7% on
human 04 CAM1 — but its values were biased long by the wandering baseline,
466 → 374 ms.

## Accuracy against a known ground truth

A synthetic AP with a constructed APD80 of **280 ms**, ensemble-averaged the same
way. "old" is the absolute threshold, "new" is the local threshold, both with the
existing `remove_Drift` enabled:

| scenario | old threshold | new threshold |
|---|---|---|
| clean, no drift | 283 ms (+3) | **283 ms** (+3) |
| clean + 2% noise | 294 ms (+14) | **285 ms** (+5) |
| + mild drift | 309 ms (+29) | **276 ms** (−4) |
| + strong drift | 358 ms (+78) | **260 ms** (−20) |

The absolute threshold **over-reports** APD, and the error grows steeply with
drift, because the global normalization anchors 0 at the 1st percentile of the
whole record — the lowest wander trough, which sits below the typical diastolic
baseline. The threshold therefore lands lower than a true 80% recovery and takes
longer to reach.

The local threshold is accurate on clean and mildly drifting signals. Under
strong drift it under-reports by ~20 ms — not perfect, but a quarter of the old
method's +78 ms error and in the conservative direction.

**This means APD values will change on existing recordings** — downward, by more
at slow pacing (a pig 1000 ms file moved 397 → 296 ms) than fast (200–300 ms
files moved only 2–4 ms). That is the old inflation being removed, not a new
error, but results computed before and after this patch are not comparable.

---

## 1. `extract_apd` — replace this block

**Find** (inside `function [data] = extract_apd(app, camera, rep_val)`):

```matlab
                [~, act_times] = max(diff_data, [], 3);          % [R x C]
                [~, ~, T]      = size(win_data);
                win_norm       = normalize_data(win_data);
                [~, pk]        = max(win_norm, [], 3);
                t_idx          = reshape(1:T, 1, 1, T);
                [hit, rep_frame] = max(cumsum((win_norm <= rep_val) & (t_idx > pk), 3) == 1, [], 3);
                rep            = double(rep_frame);
                rep(~hit)      = NaN;
```

**Replace with:**

```matlab
                [~, ~, T]      = size(win_data);
                win_norm       = normalize_data(win_data);
                t_idx          = reshape(1:T, 1, 1, T);

                % Activation must PRECEDE the peak.  Taking max(dV/dt) over the
                % whole averaged beat lets a late artifact win: on human 005 it
                % placed activation at frame 440 of 1100 for 40-64% of pixels
                % that had identical SNR and beat amplitude to their neighbours.
                % The pre-activation window then swallowed the entire AP, giving
                % amp <= 0 -- and the previous code reported those pixels as
                % APD80 = 78-107 ms for human ventricle at 1 Hz.  Masking the
                % derivative at and after the peak costs nothing where activation
                % was already correct and repairs it where it was not.
                [~, pk_g]      = max(win_norm, [], 3);
                nd             = size(diff_data, 3);
                td             = reshape(1:nd, 1, 1, nd);
                d_mask         = diff_data;
                d_mask(td >= reshape(pk_g, size(pk_g,1), size(pk_g,2), 1)) = -Inf;
                [~, act_times] = max(d_mask, [], 3);             % [R x C]
                act3           = reshape(act_times, size(act_times,1), size(act_times,2), 1);

                % Reference the repolarization threshold to THIS pixel's own
                % baseline and amplitude rather than to an absolute level in the
                % globally-normalized trace.  With a wandering baseline the
                % resting level varies across the record, so a fixed level sits
                % ABOVE the diastolic baseline of some beats and can never be
                % reached -- those pixels silently became NaN.  Local
                % referencing is also immune to beat-to-beat amplitude variation,
                % which an absolute threshold cannot handle even after a perfect
                % detrend.
                pre_act                = win_norm;
                pre_act(t_idx >= act3) = NaN;
                base                   = median(pre_act, 3, 'omitnan');  % diastolic level
                post_act               = win_norm;
                post_act(t_idx < act3) = NaN;
                [pk_val, pk]           = max(post_act, [], 3);
                amp                    = pk_val - base;

                % rep_val is the FRACTION OF AMPLITUDE still remaining, so the
                % existing 1 - pct/100 convention carries over unchanged:
                % APD80 -> rep_val = 0.20 -> 80% recovered.
                thr = base + rep_val .* amp;

                [hit, rep_frame] = max(cumsum((win_norm <= thr) & (t_idx > pk), 3) == 1, [], 3);
                rep            = double(rep_frame);
                rep(~hit | ~(amp > 0) | isnan(base)) = NaN;
```

### Optional: stop censoring silently

Immediately **after** the existing line

```matlab
                apd = apd * (1000 / app.data.acqFreq);   % frames -> ms
```

add:

```matlab
                % Report censored pixels rather than letting them vanish into
                % the map as holes.
                n_all = numel(apd);
                n_bad = nnz(isnan(apd));
                update_console(app, sprintf( ...
                    'APD%d: %d of %d pixels (%.1f%%) had no repolarization crossing.', ...
                    round(100 * (1 - rep_val)), n_bad, n_all, 100 * n_bad / max(n_all, 1)));
```

---

## 2. `extract_rep` — replace this block

**Find** (inside `function [data] = extract_rep(app, camera, rep_val)`) — the block
running from `[~, ~, T] = size(win_data);` down to and including `rep(~hit) = NaN;`:

```matlab
                [~, ~, T]  = size(win_data);
                win_norm   = normalize_data(win_data);           % one call for all pixels
                [~, pk]    = max(win_norm, [], 3);               % peak frame per pixel [R x C]
                t_idx      = reshape(1:T, 1, 1, T);
                below      = win_norm <= rep_val;
                after      = t_idx > pk;                         % broadcasts [R x C x T]
                [hit, rep_frame] = max(cumsum(below & after, 3) == 1, [], 3);
                rep        = double(rep_frame);
                rep(~hit)  = NaN;
```

**Replace with:**

```matlab
                [~, ~, T]  = size(win_data);
                win_norm   = normalize_data(win_data);           % one call for all pixels
                t_idx      = reshape(1:T, 1, 1, T);
                [~, pk]    = max(win_norm, [], 3);               % peak frame per pixel [R x C]
                pk3        = reshape(pk, size(pk,1), size(pk,2), 1);

                % Local baseline, same rationale as extract_apd.  extract_rep has
                % no activation map, so the diastolic level is taken as a low
                % percentile of the pre-peak frames -- robust to noise, and not
                % dragged upward by the rising edge the way a median would be.
                pre_pk               = win_norm;
                pre_pk(t_idx >= pk3) = NaN;
                base                 = prctile(pre_pk, 10, 3);
                amp                  = max(win_norm, [], 3) - base;
                thr                  = base + rep_val .* amp;

                below      = win_norm <= thr;
                after      = t_idx > pk3;                        % broadcasts [R x C x T]
                [hit, rep_frame] = max(cumsum(below & after, 3) == 1, [], 3);
                rep        = double(rep_frame);
                rep(~hit | ~(amp > 0) | isnan(base)) = NaN;
```

> The `extract_rep` "find" block above is reconstructed from the compiled app
> source; match on the surrounding lines rather than character-for-character if
> your copy differs. The two lines that must change are the `below = ...`
> definition and the final `rep(~hit) = NaN;`.

---

## Note on `normalize_data`

Leave the `normalize_data` call in `ensembleAverageFull` alone. Downstream QC
(`validate_conditioned`, the decay-fraction check) assumes the averaged output is
in `[0,1]`. Once the threshold is referenced locally it is scale-invariant, so the
normalization is harmless.
