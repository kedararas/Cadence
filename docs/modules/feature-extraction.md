# Feature Extraction

Extracts electrophysiological features from conditioned recordings.

- **Input:** conditioned `cmos_all_data.mat`.
- **Output:** `*-metrics.mat` (camera stacks + `ep_metrics` with all features).
- **Next step:** [Signal Analysis](signal-analysis.md) /
  [Conduction Velocity](conduction-velocity.md).

!!! tip "Video"
    TODO: embed the Feature Extraction tutorial.

## Workflow
1. **Select the input directory** of conditioned `.mat` files and **load** a file.
2. **Right-click** the image to display a representative signal; drag to find a
   good pixel. Click **Select Window** to draw a signal window over ~one cycle
   (you can add multiple windows).
3. **Select the features** to extract from the drop-downs.
4. **Select the output directory** and **Extract Features** (figures are generated).
5. **Save** the results to `*-metrics.mat`.

Batch processing applies the current feature selections to every file in the
directory.

## Features

### Voltage (Vm)
| Feature | Notes |
|---|---|
| Activation time (LAT) | At 50% upstroke |
| Repolarization time | At a chosen % |
| APD | Action-potential duration at a chosen % |
| Rise time | Between two amplitude thresholds |
| Conduction velocity | Local CV from the activation map |
| Voltage–calcium delay | EC-coupling latency (needs a Vm and a Ca camera) |

### Calcium (Ca)
Transient duration, decay time, rise time, decay τ.

!!! warning "Calcium τ at fast rates"
    See the [Signal Conditioning](signal-conditioning.md) note — measure calcium
    decay/τ from beat-windowed data, not the ensemble average, when transients
    don't relax within a cycle.

### Alternans & arrhythmia dynamics
APD/calcium alternans; spectral complexity (DF/RI/OI); wavefront tracking;
phase-singularity (rotor) dynamics.

## Notes
- TODO: recommended % thresholds per feature; pixel-selection tips; FOV entry.
