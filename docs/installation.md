# Installation

CADENCE is distributed as a standalone application — **no MATLAB license is
required** to run it.

## System requirements

| | Minimum | Recommended | Heavy (4-camera / batch) |
|---|---|---|---|
| RAM | 16 GB | 32 GB | 64 GB |
| CPU | 4-core 64-bit | 8-core | 8–16 core |
| Storage | SSD, 50 GB free | SSD, 250 GB+ | NVMe SSD, 500 GB+ |
| Display | 1920×1080 | 2560×1440 | 2560×1440+ |
| GPU | not required | not required | not required |

RAM is the main constraint (recordings are large and conditioning holds several
in-memory copies). An SSD is strongly recommended.

## macOS

1. Open `CADENCE_Installer.dmg` and run the installer inside.
2. It installs CADENCE and the **MATLAB Runtime** (free, no license).
3. Launch **CADENCE**.

!!! warning "Gatekeeper"
    If the app is not yet notarized, macOS may show "unidentified developer."
    Right-click the app → **Open** → **Open** (once). TODO: remove this note once
    the build is notarized.

## Windows

TODO: add once the Windows build is available. Run `CADENCE_Installer.exe`; it
installs the app + the Windows MATLAB Runtime.

## Running from source (developers)

Requires MATLAB + Signal Processing, Image Processing, Statistics & Machine
Learning, Curve Fitting, and Computer Vision toolboxes. Add the repo to the path
and run `Cadence.mlapp`.
