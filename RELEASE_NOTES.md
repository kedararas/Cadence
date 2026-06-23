# Release notes

Newest release first. Copy the latest entry into the GitHub Release description
when you publish, and attach the platform installer/DMG as a release asset.

---

## v2026.1 — macOS (Apple Silicon)

CADENCE is an end-to-end toolkit for processing and analyzing cardiac
optical-mapping recordings (transmembrane voltage and intracellular calcium):
data conversion, signal conditioning, feature extraction (activation, APD,
conduction velocity, calcium kinetics, alternans, arrhythmia/rotor dynamics),
visualization, and analysis. No MATLAB license required.

> **`CADENCE.dmg` is built for macOS Apple Silicon (M1/M2/M3/M4) only.** It will not run on Intel Macs.

**Download & install (macOS, Apple Silicon)**
1. Install the free **MATLAB Runtime R2025b for macOS (Apple silicon)** — on the download page, select release **R2025b** and the **macOS Apple silicon** variant (not Intel):
   https://www.mathworks.com/products/compiler/matlab-runtime.html
2. Download **`CADENCE.dmg`** (attached below), open it, and drag **CADENCE** to Applications.
3. Launch CADENCE.

The macOS build is **code-signed and notarized by Apple** — it opens normally,
no Gatekeeper warning.

**Requirements:** macOS (Apple Silicon); 16 GB RAM minimum (32 GB recommended);
SSD strongly recommended. Full system requirements and usage:
https://kedararas.github.io/Cadence/

**Notes**
- Windows and Linux builds are planned and will be released separately.
- Recordings are large; ample RAM and an SSD are recommended.

**License:** MIT — © 2026 Kedar Aras (ARAS Lab). Bundled third-party components
and the MATLAB Runtime retain their own licenses (see the repository).
