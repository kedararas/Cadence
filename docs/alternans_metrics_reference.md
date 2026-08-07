# CADENCE Alternans Metrics Reference

Every alternans and substrate metric produced by the 3-D (spatial-map) analysis
pipeline, with its output field, units, definition, and physiological
significance. Field names are the struct fields of each function's 3-D output.

**Sources:** `analyzeAPAlternans.m` · `analyzeCaTransientAlternans.m` · `assess_arrhythmia_substrate.m`

## Conventions

- All maps are `[R × C]` pixel grids, masked to tissue (NaN outside).
- **ALT** (alternans magnitude) = `mean(odd − even)` over beat pairs.
- **Ratio** = `|ALT| / mean(metric)` — dimensionless fraction.
- **Phase**: `+1` odd beats larger/longer, `−1` even beats larger/longer, `0` invalid.
- Significance defaults: alternans **ratio > 0.10** = present, **> 0.20** = strong;
  per-pixel **p < 0.05**. The composite risk map saturates the alternans component at 0.20.
- Ratio maps of the same kind (e.g. APD vs CaT amplitude ratio) are dimensionless and
  may share a colour scale; duration and amplitude maps must not.

---

## Voltage — Action-Potential Alternans

`analyzeAPAlternans.m` (3-D mode)

| Metric | Field(s) | Units | Definition | Significance | Threshold / flag |
|---|---|---|---|---|---|
| APD alternans | `APD30_alt`, `APD50_alt`, `APD80_alt` | ms | Beat-pair difference in APD at 30/50/80 % repolarisation. | Canonical repolarisation alternans and primary dispersion substrate. APD80 standard; auto-falls back to APD50/30 when unreachable at fast pacing. | ratio > **0.10** |
| APD alternans ratio | `APD50_ratio_map` | — | \|APD alternans\| / mean APD per pixel. | Magnitude normalised to local APD; feeds composite risk (saturates at 0.20). | >0.10 present · >0.20 strong |
| AP amplitude alternans | `amp_alt_map`, `amp_ratio_map` | a.u. / — | Beat-pair difference in AP amplitude (V<sub>max</sub> − V<sub>rest</sub>). | Full-tissue reference map. Less specific — sensitive to signal / mechanical artefact. | interpret with SNR |
| Upstroke (dV/dt) alternans | `dvdt_alt_map`, `dvdt_ratio_map` | a.u./ms | Beat-pair difference in max upstroke velocity. | Excitability / conduction alternans; electrical correlate of the conduction panel; driver of 2:1 block. | — |
| Triangulation alternans | `tri_alt_map` | ms | Beat-pair difference in triangulation (APD90 − APD30). | Repolarisation-morphology instability; proarrhythmic (reduced repolarisation reserve). | — |
| Alternans phase | `phase_map` | +1 / −1 / 0 | Sign of APD alternans (+1 odd longer, −1 even longer). | Building block of concordance/discordance; anti-phase regions = discordant alternans. | [−1, 1] scale |
| Continuous phase angle | `phase_angle_map` | rad (−π…π) | Angle of FFT coeff at 0.5 cyc/beat. | Spatial gradients reveal travelling alternans waves; ~π apart = anti-phase. | — |
| Spectral alternans index | `spectral_map` | power | Power at 0.5 cyc/beat in per-beat APD series. | Frequency-domain magnitude, robust to missed/extra beats. Needs ≥4 beats. | — |
| Statistical significance | `pval_map`, `tstat_map` | p / t | Pixel-wise paired t-test, odd vs even APD. | Separates sustained alternans from noise; gate low-SNR pixels. Needs ≥3 pairs. | **p < 0.05** |
| AP load / release alternans | `AP_release_alt_map`, `AP_load_alt_map`, `AP_amp_phase_map`, `AP_L_vamp_map`, `AP_S_vamp_map` | — / a.u. | Wang 2014 applied to AP amplitude: release = 1 − S/L; load = D/L; L/S = large/small-beat amplitude. | Parallels the calcium load/release decomposition for AP-vs-Ca mechanism comparison. | release ∈ [0,1] |

---

## Calcium — Ca²⁺-Transient Alternans

`analyzeCaTransientAlternans.m` (3-D mode)

| Metric | Field(s) | Units | Definition | Significance | Threshold / flag |
|---|---|---|---|---|---|
| CaT amplitude alternans | `amp_alt_map`, `amp_ratio_map` | a.u. / — | Beat-pair difference in Ca-transient amplitude (peak − onset), and its ratio. | Canonical calcium alternans and most-reported metric — SR Ca²⁺ cycling instability. Pairs with APD alternans. | ratio > **0.10** |
| CaT duration alternans | `D50_alt_map`, `D80_alt_map`, `D50_ratio_map`, `D80_ratio_map` | ms / — | Beat-pair difference in CaTD to 50 % / 80 % decay. | Calcium analogue of APD alternans. **Reachability caveat:** at fast pacing the transient may not decay to 80 % within the beat window — check `level_coverage` before trusting D80. | see `level_coverage` |
| Time-to-peak alternans | `TTP_alt_map`, `TTP_ratio_map` | ms / — | Beat-pair difference in onset-to-peak time. | Release-kinetics alternans — variability in SR Ca²⁺ release speed. | — |
| SR Ca release alternans | `release_alt_map` | — [0,1] | Wang 2014: 1 − S/L (L/S = mean large/small-beat amplitude, per-pixel). | Fractional-release instability — the release-driven component. Per-pixel L/S keeps it valid under discordant alternans. | 0 = none · →1 strong |
| SR Ca load alternans | `load_alt_map` | — (≥0) | Wang 2014: D/L, D = mean \|diastolic<sub>L</sub> − diastolic<sub>S</sub>\|. | Diastolic-load instability. **load ≈ 0 with release > 0 ⇒ release-driven** alternans (common arrhythmogenic route). | ≈0 ⇒ release-driven |
| Large / small beat amplitude | `L_mean_map`, `S_mean_map` | a.u. | Per-pixel mean amplitude of locally larger (L) / smaller (S) beats. | Underlie release/load ratios; L ≥ S everywhere by construction. | — |
| Ca alternans phase | `phase_map`, `phase_angle_map` | +1/−1/0 · rad | Sign and continuous angle of Ca amplitude alternans. | Multiply by AP `phase_map` for Ca–Vm concordance: product +1 in-phase, −1 out-of-phase. | [−1, 1] scale |
| Spectral & significance | `spectral_map`, `pval_map`, `tstat_map` | power · p/t | Power at 0.5 cyc/beat; paired t-test odd vs even amplitude. | Same roles as the voltage versions. | **p < 0.05** |
| CaTD reachability | `level_coverage` | fraction | Fraction of tissue-beats with a finite CaTD per decay level (diagnostic). | Data-quality gate: <0.5 at D80 means that map is mostly NaN — don't over-interpret. | warn if **< 0.50** |

---

## Substrate & Ca–AP Coupling

`assess_arrhythmia_substrate.m`

| Metric | Field(s) | Units | Definition | Significance | Threshold / flag |
|---|---|---|---|---|---|
| Concordance / discordance | `concordance_ratio`, `is_discordant` | fraction · bool | Fraction of significant pixels sharing the majority phase. | Discordant alternans steepens repolarisation gradients — classic precursor to block and reentry. | discordant if **< 0.80** |
| Nodal lines | `nodal_lines` | logical | Pixels bordering both +phase and −phase regions. | Phase-reversal boundaries — highest-risk sites where unidirectional block initiates. | binary mask |
| Repolarisation dispersion | `apd_disp_global`, `gradient_mag`, `max_gradient` | ms · ms/px | SD of mean-APD map and its Sobel gradient magnitude/peak. | Large dispersion / steep gradients mark repolarisation borders — block lines and reentry substrate. | — |
| Ca–AP coupling | `coupling_map`, `coupling_fraction`, `inphase_mask` | r ∈ [−1,1] | Per-pixel Pearson r between per-beat Ca-amplitude and APD series. | r > 0 in-phase (Ca drives AP via NCX / I<sub>CaL</sub>); r < 0 out-of-phase (voltage drives Ca). Out-of-phase = higher-risk substrate. | >0.6 Ca-driven · <0.4 V-driven |
| Ca metrics (comparison) | `ca_phase_map`, `ca_phase_angle_map`, `ca_alt_ratio_map` | +1/−1/0 · rad · — | Calcium phase and alternans-ratio maps wired into the substrate struct. | Overlay Ca and AP phase/magnitude on the same significant-pixel footprint. | — |
| Diastolic Ca elevation | `diast_ca_map`, `diast_ca_alt_map` | norm. a.u. | Mean normalised diastolic Ca (last 10 % of interval) and its odd-even alternans. | Elevated diastolic Ca → DADs / triggered activity (independent risk axis); its alternans indicates SR-loading instability. | requires CaData |
| APD restitution slope | `restitution_map`, `slope_gt1_mask`, `mean_slope` | dAPD/dDI | Local slope of APD vs diastolic interval from odd/even (APD, DI) pairs. | Restitution-hypothesis instability predictor: slope > 1 dynamically unstable, alternans/VF-susceptible. | unstable if **> 1** |
| Composite arrhythmia risk | `risk_map`, `risk_global`, `risk_components` | [0, 1] | Weighted sum (norm. to 1): alternans ratio (1.0) + APD gradient (1.0) + nodal proximity (1.0) + restitution>1 (1.0) + Ca out-of-phase (0.5) + diastolic Ca (0.5). | Single per-pixel susceptibility score integrating dispersion, dynamic instability, and Ca-driven amplification. | higher = higher risk |

---

### Notes for reporting

- **AP amplitude units:** the 3-D maps operate on optical fluorescence (F/F₀-style), so AP
  amplitude is reported in arbitrary units (a.u.), not mV, even though the 1-D helper docs mention mV.
- **Significance cut-offs** (ratio > 0.10, p < 0.05) are the pipeline defaults; adjust if your
  analysis runs used different thresholds.
