#!/bin/bash
# Sign + notarize + staple BOTH macOS builds (Apple Silicon and Intel) in sequence.
# Run from anywhere; it cd's to the project root. Each leg waits on Apple, so the
# whole thing takes several minutes. Edit the paths below if your build folders differ.

set -euo pipefail
cd "$(dirname "$0")/.."   # project root

echo "==================  APPLE SILICON  =================="
./notarize/sign_notarize.sh \
  distribution/macos-apple-silicon/CADENCE_AppleSilicon.app \
  distribution/macos-apple-silicon/CADENCE-AppleSilicon.dmg  CADENCE

echo
echo "==================  INTEL  =================="
./notarize/sign_notarize.sh \
  distribution/macos-intel/CADENCE_Intel.app \
  distribution/macos-intel/CADENCE-Intel.dmg  CADENCE

echo
echo ">> Both done. Notarized distributables:"
echo "   distribution/macos-apple-silicon/CADENCE-AppleSilicon.dmg"
echo "   distribution/macos-intel/CADENCE-Intel.dmg"
