# Conduction Velocity

Computes local and directional conduction velocity from activation maps.

- **Input:** `*-metrics.mat` (from Feature Extraction).
- **Output:** CV values and maps.

!!! tip "Video"
    TODO: embed the Conduction Velocity tutorial.

## Workflow
1. **Select the input directory** of `*-metrics.mat` files and load a recording.
2. **Enter the field-of-view (FOV) information.** Pixel resolution = camera
   sensor size ÷ pixel count (e.g. SciMedia N256: 17.6 mm ÷ 256 ≈ 0.07 mm/px).
3. **Select the points of interest** — typically the pacing site (Point A) and
   the location where CV is measured (Point B).
   - CV is most accurate within ~20 mm of the pacing site, with ≥1 mm/px resolution.
4. **Select the CV method:** Euclidean distance, Single vector, Average (multi)
   vector, or All of the above (for comparison).
5. **Calculate Conduction Velocity.**

## Methods
| Method | Idea |
|---|---|
| Euclidean (two-point) | distance/time between origin and a measurement point, along/across the fiber axis |
| Multi-vector (directional) | pools the local CV vector field over all valid tissue; reports longitudinal vs transverse CV |

!!! note "Planar wavefronts"
    For a uni-directional (planar) sweep, longitudinal vs transverse CV is not
    meaningful — CADENCE warns when wavefront planarity is high. Trust only the
    along-propagation value in that case.

## Notes
- TODO: bin-size selection and its effect on the local-CV kernel.
- TODO: interpreting the CV map / statistics.
