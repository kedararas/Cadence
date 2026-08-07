# Troubleshooting

!!! note
    Each module's console (the text panel) reports what it's doing and any
    errors. Check it first.

## Loading / data
- **"SCHEMA FAIL …" on load** — the `.mat` is missing a required field or has the
  wrong type for the stage. Re-run the producing module. See [Data format](data-format.md).
- **Recording skipped in a batch** — usually a schema FAIL on that file; the batch
  continues with the rest.

## Conditioning
- **Ensemble average looks wrong (doesn't decay)** — expected when the signal
  doesn't relax within one pacing cycle (e.g. calcium at fast rates). Use
  beat-windowed data for decay/τ features, not the average.

## Feature extraction
- **A feature "could NOT be extracted"** — TODO: common causes (no signal window,
  pixel on non-tissue, threshold too strict).

## Conduction velocity
- **CV warning about a planar / uni-directional wavefront** — longitudinal vs
  transverse CV isn't meaningful for a planar sweep; trust the along-propagation value.

## Installation / launch
- **macOS "unidentified developer"** — right-click the app → Open → Open (once),
  until the build is notarized.
- **App won't launch** — confirm the MATLAB Runtime installed by the installer
  matches the build version.

---

TODO: add FAQ entries as users report issues. Link a GitHub Issues template here.
