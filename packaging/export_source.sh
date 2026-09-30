#!/bin/bash
# Export a clean, portable CADENCE folder for INTERNAL distribution: MATLAB
# source + app images (including the internal batch tools), and the released
# installers. Nothing else from the working directory is copied. The public
# gets the GitHub release instead.
#
# Usage (from anywhere):
#   packaging/export_source.sh                       # HEAD -> ~/Desktop/CADENCE_v<version>
#   packaging/export_source.sh --ref v1.0.1          # source exactly as tagged
#   packaging/export_source.sh --out /Volumes/USB    # write somewhere else
#   packaging/export_source.sh --no-batch            # leave out the internal batch tools
#   packaging/export_source.sh --no-installers       # source only
#
# Output
#   CADENCE_v<version>/
#     source/       the MATLAB root: the six .mlapp apps, every helper folder they
#                   use, the internal batch tools (batch_extraction_helper), the
#                   validation harness (manual_validation_helper, cited in the
#                   paper), the images the apps and builds use, packaging/
#                   build scripts, and the license / citation / README files.
#                   Open MATLAB here and run:  addpath(genpath(pwd)); Cadence
#     installers/   CADENCE-<version>-AppleSilicon.dmg and CADENCE_Installer_Windows.exe
#                   from distribution/, if present
#     EXPORT_INFO.txt  version, git commit, date, and contents
#   CADENCE_v<version>.zip   the same folder, zipped
#
# Only files committed to git are exported (git archive), so build output
# (distribution/, release/), Claude worktrees (.claude/), the MATLAB project
# metadata (resources/, *.prj), Illustrator sources (*.ai), docs and signing
# scripts never reach the copy. Uncommitted edits are NOT included; commit first.

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
REF="HEAD"
OUT="$HOME/Desktop"
WITH_BATCH=1
WITH_INSTALLERS=1
while [ $# -gt 0 ]; do
  case "$1" in
    --ref) REF="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --no-batch) WITH_BATCH=0; shift ;;
    --with-batch) WITH_BATCH=1; shift ;;   # default; kept for compatibility
    --no-installers) WITH_INSTALLERS=0; shift ;;
    -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

cd "$REPO"
COMMIT="$(git rev-parse --short "$REF")"
VERSION="$(git show "$REF:CITATION.cff" | sed -n 's/^version: *"\([^"]*\)".*/\1/p')"
[ -n "$VERSION" ] || { echo "Could not read version from CITATION.cff at $REF" >&2; exit 1; }
if [ "$REF" = "HEAD" ] && [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "NOTE: uncommitted changes exist; they are NOT in the export (it uses committed files)." >&2
fi

NAME="CADENCE_v$VERSION"
DEST="$OUT/$NAME"
[ -e "$DEST" ] && { echo "ERROR: $DEST already exists; remove or rename it first." >&2; exit 1; }

# ---- what goes into source/ -------------------------------------------------
PATHS=(
  Cadence.mlapp
  Cadence_Data_Conversion.mlapp
  Cadence_Signal_Conditioning.mlapp
  Cadence_Feature_Extraction.mlapp
  Cadence_Signal_Analysis.mlapp
  Cadence_Conduction_Velocity.mlapp
  utils
  signal_conditioning_helper
  feature_extraction_helper
  conduction_velocity_helper
  arrhythmia_dynamics_helper
  validation_helper
  manual_validation_helper
  Logo_v3.png Sidebar.png Sidebar_v1.jpg assets
  packaging
  License.txt THIRD_PARTY_LICENSES.md README.md RELEASE_NOTES.md CITATION.cff
  WINDOWS_BUILD_CHECKLIST.md
)
[ "$WITH_BATCH" = 1 ] && PATHS+=(batch_extraction_helper)

mkdir -p "$DEST/source"
git archive --format=tar "$REF" -- "${PATHS[@]}" | tar -x -C "$DEST/source"

# ---- installers --------------------------------------------------------------
INST=()
if [ "$WITH_INSTALLERS" = 1 ]; then
  mkdir -p "$DEST/installers"
  for f in "distribution/macos-apple-silicon/CADENCE-$VERSION-AppleSilicon.dmg" \
           "distribution/windows/CADENCE_Installer_Windows.exe"; do
    if [ -f "$f" ]; then cp "$f" "$DEST/installers/"; INST+=("$(basename "$f")")
    else echo "WARNING: installer not found, skipped: $f" >&2; fi
  done
fi

# ---- record what this is -------------------------------------------------------
{
  echo "CADENCE v$VERSION — portable export"
  echo "Exported : $(date '+%Y-%m-%d %H:%M')"
  echo "Git ref  : $REF (commit $COMMIT)"
  echo "Batch tools included: $([ "$WITH_BATCH" = 1 ] && echo yes || echo no)"
  echo
  echo "source/      MATLAB root. In MATLAB: cd to this folder, then"
  echo "               addpath(genpath(pwd)); Cadence"
  echo "             Requires MATLAB R2025b + Signal Processing, Image Processing,"
  echo "             Statistics and Machine Learning, Curve Fitting toolboxes."
  echo "             Windows build: cd packaging; build_windows('$VERSION')"
  echo "installers/  ${INST[*]:-(none)}"
  echo "             macOS: install MATLAB Runtime R2025b (Apple silicon), open the DMG."
  echo "             Windows: run the installer (downloads the Runtime; unsigned ->"
  echo "             SmartScreen: More info -> Run anyway)."
  echo
  echo "Files in source/:"
  (cd "$DEST/source" && find . -type f | sed 's|^\./|  |' | sort)
} > "$DEST/EXPORT_INFO.txt"

(cd "$OUT" && ditto -c -k --keepParent "$NAME" "$NAME.zip")

echo "Exported $NAME from $REF ($COMMIT)"
echo "  folder: $DEST"
echo "  zip   : $OUT/$NAME.zip ($(du -h "$OUT/$NAME.zip" | cut -f1))"
echo "  source files: $(find "$DEST/source" -type f | wc -l | tr -d ' ')   installers: ${#INST[@]}"
