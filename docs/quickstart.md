# Quick start

A first end-to-end run, from raw data to results. Each step links to the full
module page.

!!! tip "Video"
    TODO: embed the "CADENCE in 10 minutes" overview tutorial here.

## 1. Convert raw data → `.mat`
Open **[Data Conversion](modules/data-conversion.md)**, select the input
directory of raw recordings, set the sampling frequency, choose an output
directory, and convert.

## 2. Condition the signals
Open **[Signal Conditioning](modules/signal-conditioning.md)**, select the
converted `.mat` files, choose conditioning options, and process (single or
batch). Produces conditioned `.mat`.

## 3. Extract features
Open **[Feature Extraction](modules/feature-extraction.md)**, load a conditioned
file, pick a representative pixel and signal window, select the features to
extract, run, and save `*-metrics.mat`.

## 4. Analyze & visualize
Open **[Signal Analysis](modules/signal-analysis.md)** or
**[Conduction Velocity](modules/conduction-velocity.md)** and load the
`*-metrics.mat` file to visualize, curate, quantify, and export.

---

TODO: add a screenshot of the launcher and a one-paragraph "what you'll need"
(a sample dataset link, expected runtime).
