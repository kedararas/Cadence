# CADENCE — Windows Build & Signing Checklist

Building CADENCE for Windows is a **separate build on a Windows machine** — the
macOS app/installer are Apple-Silicon binaries and do not run on Windows, and
MATLAB Compiler does not cross-compile. There is **no notarization on Windows**
(that's Apple-only); the equivalent is **Authenticode code signing**.

## 1. Prerequisites (on a Windows PC)
- Windows 10/11, 64-bit.
- **MATLAB R2025b** (match the macOS build version) + **MATLAB Compiler**.
- Toolboxes: **Signal Processing, Image Processing, Statistics & Machine Learning,
  Curve Fitting, Computer Vision**.
- For signing: a Windows **Authenticode code-signing certificate** (+ its hardware
  token / HSM) and **`signtool.exe`** (from the Windows SDK).

## 2. Copy the project to Windows
Bring the whole project, same folder structure: `Cadence.mlapp` (launcher), the
five module `.mlapp` files, all helper folders (`utils/`, `*_helper/`,
`validation_helper/`), `Logo_v3.png`, `Sidebar.png`, the icon, and the docs.

## 3. Build (Application Compiler)
- Run `applicationCompiler`.
- **Main file:** `Cadence.mlapp`.
- **Files installed with the app:** `Logo_v3.png`, `Sidebar.png`, icon.
- **Application name:** `CADENCE`; **Version:** match the macOS build.
- **Application icon:** a `.ico` (convert the icon PNG to `.ico`).
- **Additional Installer Options → Installer name:** `CADENCE_Installer`.
- **Runtime:** web (downloads at install) or included (offline, larger).
- **Save the deployment project** (e.g. `CadenceDeployWin.prj`) so rebuilds are one click.
- **Package** → produces `CADENCE.exe` + a Windows installer (`CADENCE_Installer` setup `.exe`).
- Verify a clean build: `unresolvedSymbols.txt` empty, nothing excluded (same checks as macOS).

## 4. Code signing — DECIDE THE CERTIFICATE FIRST (lead time + cost)
Signing is recommended: an unsigned `.exe`/installer triggers a **SmartScreen**
"Windows protected your PC" warning (and some antivirus flags). It still runs, but
looks alarming.

**Certificate decision:**
- **OV (Organization Validation)** — cheaper (~$100–300/yr). SmartScreen reputation
  builds up over downloads/time, so early users may still see a warning.
- **EV (Extended Validation)** — **immediate** SmartScreen trust (no warning from
  the first download); pricier (~$300–600/yr).
- All publicly-trusted code-signing certs now must be stored on a **hardware token
  / HSM** (the CA ships one, or use a cloud HSM) — you can't keep a `.pfx` on disk.
- CAs: DigiCert, Sectigo, GlobalSign, SSL.com, etc.
- The cert is issued to an **organization** (it gets validated). This ties to the
  same "who owns CADENCE" question as the Apple account — University at Buffalo as
  the cert holder vs an individual/sole-proprietor cert. Sort this before purchase.

**Sign both the app and the installer** (Windows SDK `signtool`), with an RFC-3161
timestamp so signatures stay valid after the cert expires:
```bat
signtool sign /fd SHA256 /tr http://timestamp.digicert.com /td SHA256 /a CADENCE.exe
signtool sign /fd SHA256 /tr http://timestamp.digicert.com /td SHA256 /a CADENCE_Installer.exe
signtool verify /pa /v CADENCE_Installer.exe
```
(`/a` auto-selects the cert; for an HSM use the token/CSP options the CA provides.)

## 5. Test on a clean Windows machine
- A Windows box with **no MATLAB**. Run the installer (installs the app + the
  **Windows** MATLAB Runtime R2025b — a separate download from the macOS Runtime).
- Launch CADENCE, run the pipeline end-to-end with a real `.mat`.
- Confirm SmartScreen behavior matches the cert type (no warning for EV).

## 6. Per release
Re-sign each build you distribute (signing is per-binary, like macOS). **No
notarization step.** Windows distributes the signed installer directly (no DMG).

## Quick macOS vs Windows reference
| | macOS | Windows |
|---|---|---|
| Build | App Compiler → `.app` | App Compiler → `.exe` |
| Sign | `codesign` (Apple Developer ID) | `signtool` (Authenticode cert) |
| Notarize | Yes (Apple) | None |
| Gatekeeper-equivalent | Gatekeeper | SmartScreen |
| Package | `.dmg` | installer `.exe` |
| Runtime | macOS MATLAB Runtime | Windows MATLAB Runtime |
