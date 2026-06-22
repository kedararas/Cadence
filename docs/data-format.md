# Data format

CADENCE stores everything in a single MATLAB struct per recording, saved in a
`.mat` file (v7.3 / HDF5). Two struct shapes flow through the pipeline.

## `cmos_all_data` (raw / conditioned)
Written by Data Conversion and Signal Conditioning.

| Field | Meaning |
|---|---|
| `CAM<n>` | image stack for camera *n* (`rows × cols × frames`) |
| `CAM<n>_image` | preview/background image |
| `acqFreq` | acquisition frame rate (Hz) |
| `analog1` | analog/pacing channel (optional) |
| `num_files` | number of cameras |
| `CAM<n>_SNR` | per-camera SNR map (after conditioning) |
| `CAM<n>_average` | ensemble-averaged representative beat (after conditioning) |

## `metrics` (post feature extraction)
Written by Feature Extraction (`*-metrics.mat`), consumed by Signal Analysis and
Conduction Velocity.

| Field | Meaning |
|---|---|
| `CAM<n>`, `CAM<n>_image` | camera stacks + previews |
| `acqFreq`, `num_files` | as above |
| `window` | selected signal window(s) |
| `data_masks` | user/SNR masks |
| `ep_metrics` | all extracted features (activation, APD, CV, calcium, alternans, arrhythmia dynamics, …) |

!!! note "Schema validation"
    On load, each module structurally validates the struct (a "schema gate") and
    reports any missing/mis-typed fields, so malformed files are caught at the
    boundary rather than failing deep inside a stage.

TODO: document the `ep_metrics` sub-fields (per feature) if users need to read
the `.mat` directly.
