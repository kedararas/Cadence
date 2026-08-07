# CADENCE

**Cardiac optical-mapping signal processing and analysis toolkit** — from raw
acquisition through conditioning, feature extraction, and analysis of
transmembrane voltage and intracellular calcium recordings.

!!! note "This site is the living user manual"
    It is versioned alongside the code, so it stays in sync with each release.
    For the scientific methods and validation, see the publication
    (TODO: add citation/DOI). For citing the software, see `CITATION.cff`.

## The pipeline

CADENCE is five modules run in sequence; each writes a `.mat` consumed by the next.

| # | Module | Does | Output |
|---|--------|------|--------|
| 1 | [Data Conversion](modules/data-conversion.md) | Import raw `.gsd/.gsh` / TIFF | `cmos_all_data.mat` (raw) |
| 2 | [Signal Conditioning](modules/signal-conditioning.md) | Denoise, filter, normalize, ensemble-average | `cmos_all_data.mat` (conditioned) |
| 3 | [Feature Extraction](modules/feature-extraction.md) | Activation, APD, CV, calcium, alternans, rotors | `*-metrics.mat` |
| 4 | [Signal Analysis](modules/signal-analysis.md) | Visualize, curate, quantify, export | curated outputs |
| 5 | [Conduction Velocity](modules/conduction-velocity.md) | Local + directional CV | CV maps |

## Get started
- [Installation](installation.md) — install the app + MATLAB Runtime
- [Quick start](quickstart.md) — your first end-to-end run
- [Data format](data-format.md) — what's inside the `.mat` files
- [Video tutorials](tutorials.md) — ARAS Lab walkthroughs

## About
Developed by **Kedar Aras** (ARAS Lab). Released under the MIT License.
Contact: aras.research.lab@gmail.com · <https://github.com/kedararas/Cadence>
