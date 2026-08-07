# Release notes

Newest release first. Copy the latest entry into the GitHub Release description
when you publish, and attach the platform installer/DMG as a release asset.

---

## v2026.1 — macOS (Apple Silicon + Intel)

CADENCE is an end-to-end toolkit for processing and analyzing cardiac
optical-mapping recordings (transmembrane voltage and intracellular calcium):
data conversion, signal conditioning, feature extraction (activation, APD,
conduction velocity, calcium kinetics, alternans, arrhythmia/rotor dynamics),
visualization, and analysis. No MATLAB license required.

Both macOS builds are **code-signed and notarized by Apple** — they open
normally, no Gatekeeper warning.

**Pick the build for your Mac, and install the matching MATLAB Runtime R2025b**
(free, no license — choose release **R2025b** on the download page):
https://www.mathworks.com/products/compiler/matlab-runtime.html

| Your Mac | Download | MATLAB Runtime variant |
|---|---|---|
| **Apple Silicon** (M1/M2/M3/M4) | `CADENCE-AppleSilicon.dmg` | R2025b — **macOS Apple silicon** |
| **Intel** | `CADENCE-Intel.dmg` | R2025b — **macOS Intel** |

*(Not sure? `CADENCE-Intel.dmg` also runs on Apple Silicon via Rosetta 2, but the
Apple Silicon build is faster on M-series Macs.)*

**Install**
1. Install the matching **MATLAB Runtime R2025b** (table above).
2. Open the DMG for your Mac, drag **CADENCE** to Applications.
3. Launch CADENCE.

**Requirements:** macOS (Apple Silicon or Intel); 16 GB RAM minimum (32 GB
recommended); SSD strongly recommended. Full system requirements and usage:
https://kedararas.github.io/Cadence/

**Notes**
- Windows and Linux builds are planned and will be released separately.
- Recordings are large; ample RAM and an SSD are recommended.

**License:** MIT — © 2026 Kedar Aras (ARAS Lab). Bundled third-party components
and the MATLAB Runtime retain their own licenses (see the repository).
