# Signal Conditioning

Cleans and prepares raw recordings for feature extraction.

- **Input:** raw `cmos_all_data.mat` (from Data Conversion), or `.gsh` folders.
- **Output:** conditioned `cmos_all_data.mat` (adds `CAM<n>_SNR`,
  `CAM<n>_average`, and the processed stacks).
- **Next step:** [Feature Extraction](feature-extraction.md).

!!! tip "Video"
    TODO: embed the Signal Conditioning tutorial.

## Workflow
1. **Select the input directory** and pick the file(s) to process.
2. **Choose conditioning options** (each toggled independently — see below).
3. **Select the output directory** and run (single file or batch).

## Conditioning stages
Applied in this order; each is optional unless noted.

| Stage | What it does |
|---|---|
| SNR mapping | Per-camera signal-to-noise map (always on; later stages depend on it) |
| Drift correction | Removes baseline wander |
| SVD denoising | Low-rank spatial denoising (falls back to binning per camera if needed) |
| Spatial binning | N×N box filter (3/5/7/9) to raise SNR |
| Temporal filtering | Savitzky–Golay low-pass at the chosen cutoff |
| Motion correction | Registers frames to a reference beat |
| Normalization | Scales each pixel's signal to a common range |
| Ensemble averaging | Builds a representative beat (`CAM<n>_average`) |
| Inversion check | Detects and flips inverted signals |
| SNR masking | Masks non-tissue / low-SNR pixels |

!!! warning "Ensemble averaging at fast pacing rates"
    The ensemble average can be unreliable when the signal does not relax within
    one cycle (e.g. calcium transients at very fast pacing). CADENCE flags this
    (decay-fraction QC). Do **not** use the ensemble average for calcium
    decay/τ/duration features on such recordings — use beat-windowed data.

## Notes
- TODO: recommended option presets per species / recording type.
- The full beat-to-beat stack (`CAM<n>`) is always preserved; ensemble averaging
  only adds `CAM<n>_average`, so alternans/arrhythmia analysis is unaffected.
