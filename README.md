# CADENCE

**C**ardiac **A**nalysis, **D**ynamics, and **EN**semble **C**haracterization **E**ngine — an end‑to‑end toolkit for processing and analyzing cardiac optical‑mapping recordings (transmembrane voltage and intracellular calcium).

> CADENCE takes raw high‑speed optical‑mapping data from acquisition through conditioning, feature extraction, and analysis — producing activation maps, action‑potential / calcium‑transient metrics, conduction velocity, alternans, and arrhythmia (wavefront / rotor) dynamics.

**Documentation:** <https://kedararas.github.io/Cadence/> &nbsp;·&nbsp; **Video tutorials:** ARAS Lab YouTube (see the docs site)

---

## Table of contents
- [Overview](#overview)
- [Documentation](#documentation)
- [Modules & workflow](#modules--workflow)
- [System requirements](#system-requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Data format](#data-format)
- [Running from source (developers)](#running-from-source-developers)
- [Citing CADENCE](#citing-cadence)
- [License](#license)
- [Third‑party components](#third-party-components)
- [Authors & contact](#authors--contact)

---

## Overview

CADENCE is a MATLAB® App Designer application distributed as a **standalone executable** — end users do **not** need a MATLAB license. It is organized as five modules that form a sequential pipeline, launched from a single CADENCE launcher. Every module validates the structure of the data it loads before processing (schema/QC gates), so malformed files are caught at the boundary rather than failing deep inside a stage.

Typical inputs are dual‑camera (voltage + calcium) or single‑camera recordings; data are large (a single recording is commonly hundreds of MB), so see [System requirements](#system-requirements).

## Documentation

The full user manual lives at **<https://kedararas.github.io/Cadence/>** — install guide, a page per module, the data format, and links to the **ARAS Lab** video tutorials. The docs are written in Markdown under [`docs/`](docs/) and **versioned with the code**, so they stay in sync with each release.

Build/preview/publish the docs site (MkDocs + Material):
```bash
pip install mkdocs-material        # one-time
mkdocs serve                       # live preview at http://127.0.0.1:8000
mkdocs gh-deploy                   # publish to GitHub Pages
```

This README is the short reference; the docs site is the living manual, and the scientific methods/validation are in the publication (see [Citing CADENCE](#citing-cadence)).

## Modules & workflow

The modules are designed to be run in order; each writes a `.mat` file consumed by the next.

```
 Raw acquisition (.gsd/.gsh, .tif)
            │
            ▼
 ┌─────────────────────┐
 │ 1. Data Conversion  │  →  cmos_all_data.mat  (raw)
 └─────────────────────┘
            │
            ▼
 ┌─────────────────────┐
 │ 2. Signal           │  →  cmos_all_data.mat  (conditioned)
 │    Conditioning     │
 └─────────────────────┘
            │
            ▼
 ┌─────────────────────┐
 │ 3. Feature          │  →  *-metrics.mat
 │    Extraction       │
 └─────────────────────┘
            │
        ┌───┴───────────────┐
        ▼                   ▼
 ┌──────────────┐   ┌────────────────────┐
 │ 4. Signal    │   │ 5. Conduction      │
 │    Analysis  │   │    Velocity        │
 └──────────────┘   └────────────────────┘
```

**1. Data Conversion** — Imports raw SciMedia `.gsd/.gsh` and TIFF stacks into the CADENCE `cmos_all_data` format: per‑camera image stacks, preview images, acquisition frame rate, and the analog/pacing channel. Supports batch processing of a directory.

**2. Signal Conditioning** — SNR mapping, baseline‑drift removal, SVD spatial denoising, spatial binning, temporal (Savitzky–Golay) filtering, motion‑artifact correction, normalization, ensemble averaging, signal‑inversion detection/correction, and SNR masking. Supports batch processing of a directory.

**3. Feature Extraction** —
- *Voltage (Vm):* local activation time (LAT), repolarization time, action‑potential duration (APD), upstroke rise time, conduction velocity, and voltage–calcium activation delay.
- *Calcium (Ca):* transient duration, decay time, rise time, and decay time constant (τ).
- *Alternans:* AP and calcium‑transient alternans.
- *Arrhythmia dynamics:* spectral complexity (dominant frequency / regularity / organization), wavefront tracking, and phase‑singularity (rotor) dynamics.
- *Image registration:* image alignment (dual Vm-Ca imaging).
- *Data masking:* SNR-based data masking, user created masks
- Supports single‑file and batch processing.

**4. Signal Analysis** — Visualize and curate feature maps; apply ROI and SNR data masks and histogram bounds; compute per‑region statistics; play phase / annotated‑rotor movies; export maps, representative signals, and movies.

**5. Conduction Velocity** — Local CV vector fields and CV estimation by the Euclidean (two‑point) and multi‑vector directional methods (longitudinal vs. transverse relative to the estimated fiber axis), with activation‑ and CV‑map display.

## System requirements

CADENCE’s requirements are driven by data size — recordings are large and the conditioning/ensemble steps hold multiple in‑memory copies, so **RAM is the dominant constraint**. A discrete GPU is **not** required or used.

### Software
- **MATLAB Runtime `R20XX?` (free, no license).** This is the only mandatory install and **must match the version CADENCE was built with**. *(Replace `R20XX?` with your build release, e.g. R2024b.)*
- **Operating system (64‑bit):** Windows 10/11, macOS (Apple Silicon and Intel are **separate builds** — install the one matching your Mac), or a Linux distribution supported by the build’s MATLAB release.
- Administrator rights for the one‑time Runtime install; ~2–4 GB free disk for the Runtime.
- **No MATLAB license and no toolboxes are required by the end user** — all dependencies are bundled.

### Hardware

| | Minimum | Recommended | Heavy (4‑camera / batch) |
|---|---|---|---|
| **RAM** | 16 GB | 32 GB | 64 GB |
| **CPU** | 4‑core 64‑bit | 8‑core, modern | 8–16 core |
| **Storage** | SSD, 50 GB free | SSD, 250 GB+ | NVMe SSD, 500 GB+ |
| **Display** | 1920×1080 | 2560×1440 | 2560×1440+ |
| **GPU** | not required | not required | not required |

An SSD is strongly recommended — loading hundreds‑of‑MB recordings from a spinning disk is slow. Requirements scale with recording length and camera count.

## Installation

CADENCE requires the free **MATLAB Runtime (R2025b)** — no MATLAB license needed. Install it once from MathWorks: <https://www.mathworks.com/products/compiler/matlab-runtime.html>

### macOS
1. Install the **MATLAB Runtime R2025b** (link above).
2. Open **`CADENCE.dmg`** and drag **CADENCE** to your Applications folder.
3. Launch **CADENCE** and begin at the Data Conversion module.

> The macOS build is **code-signed and notarized** by Apple, so it opens normally — no Gatekeeper warning.

### Windows / Linux
Coming soon — builds for these platforms are produced separately. The installer will bundle/download the matching MATLAB Runtime.

## Quick start

1. **Data Conversion** → select the input directory of raw recordings, set the sampling frequency, choose an output directory, and convert to `.mat`.
2. **Signal Conditioning** → select the converted `.mat` files, choose conditioning options, and process (single or batch). Outputs conditioned `.mat`.
3. **Feature Extraction** → load a conditioned file, pick a representative pixel and signal window, select the features to extract, run, and save `*-metrics.mat`.
4. **Signal Analysis** / **Conduction Velocity** → load `*-metrics.mat` to visualize, curate, quantify, and export.

## Data format

CADENCE stores everything in a single MATLAB struct per recording:
- **`cmos_all_data`** (raw / conditioned): `CAM<n>` image stacks, `CAM<n>_image` previews, `acqFreq`, optional `analog1` (pacing), and conditioning products (`CAM<n>_SNR`, `CAM<n>_average`).
- **`metrics`** (post‑extraction, the `*-metrics.mat` files): the camera stacks plus `ep_metrics` (all extracted features), `acqFreq`, `num_files`, signal windows, and data masks.

## Running from source (developers)

Requires **MATLAB** plus these toolboxes:
**Signal Processing**, **Image Processing**, **Statistics and Machine Learning**, **Curve Fitting**, and **Computer Vision**. *(Building the standalone additionally requires **MATLAB Compiler**.)*

1. Clone/copy the repository and add it to the MATLAB path (`addpath(genpath(pwd))`).
2. Open `Cadence.mlapp` in App Designer and **Run**, or run it from the command window.
3. To build the standalone: `applicationCompiler` (entry point `Cadence.mlapp`), or `compiler.build.standaloneApplication("Cadence.mlapp", ...)`. Build separately on each target OS (the Compiler does not cross‑compile).

## Citing CADENCE

If you use CADENCE in your research, please cite it. Citation metadata is in
[`CITATION.cff`](CITATION.cff) — on GitHub this powers the **“Cite this repository”**
button (APA/BibTeX export). Once you archive a release (e.g. on Zenodo) and obtain a
DOI, add it to `CITATION.cff` and reference it here.

## License

CADENCE is released under the **MIT License** — see [`License.txt`](License.txt). You are free to use, copy, modify, and distribute CADENCE, including in commercial and proprietary work, provided the copyright notice and permission notice are retained. The software is provided “as is”, without warranty of any kind.

```
MIT License — Copyright (c) 2026 Kedar Aras
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software ... THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
(See License.txt for the full text.)
```

## Third‑party components

CADENCE bundles a few small, permissively licensed helper functions from the MATLAB File Exchange (`real2rgb`, `rescale_sat`, `rgb`, `DiscreteFrechetDist`, `contourcs`). Each remains under its original (BSD/MIT‑style) license, which is compatible with CADENCE’s AGPL‑3.0 distribution. See [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) for the full inventory, authors, and licenses.

The bundled **MATLAB Runtime** is © The MathWorks, Inc., redistributed under MathWorks’ terms; it is a separately installed runtime required only to execute the compiled application.

## Authors & contact

**Kedar Aras** — ARAS Lab
Contact: aras.research.lab@gmail.com

Issues and contributions: https://github.com/kedararas/Cadence
