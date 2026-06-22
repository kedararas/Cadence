# Signal Analysis

Visualize, curate, quantify, and export the extracted feature maps.

- **Input:** `*-metrics.mat` (from Feature Extraction). Multiple files for one
  recording (e.g. panoramic mapping) can be loaded together.
- **Output:** curated maps, statistics, exported figures / signals / movies.

!!! tip "Video"
    TODO: embed the Signal Analysis tutorial.

## Workflow
1. **Select the input directory** of `*-metrics.mat` files and load a recording.
2. **Choose a feature** to visualize from the drop-down. Right-click the image
   to display representative signals (up to 4 per camera).
3. **Curate the maps** with data masks and histogram bounds:
   - Draw ROI masks; combine with global/regional masks.
   - Set lower/upper histogram bounds, enable **Use bounds**, and redisplay.
   - Masks and bounds are not mutually exclusive. Uncheck both and redisplay to reset.
4. **Statistical summary tab** — min/max/mean/std/median per mask (useful for
   spatial gradients: base vs apex, epi vs endo, etc.).
5. **Export** — curated map (+ histogram) as PDF, representative signals, or movies.

## Notes
- TODO: list the available feature maps and what each shows.
- TODO: movie/phase/annotated-rotor playback options.
