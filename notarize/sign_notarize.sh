#!/bin/bash
# Sign (hardened runtime + entitlements), DMG, notarize, and staple ONE app.
# Architecture-agnostic — run it once per build (Apple Silicon and Intel) from
# the project root. Notarization works for any architecture from any Mac.
#
# Usage:
#   ./notarize/sign_notarize.sh <path-to-.app> <output.dmg> [VolumeName]
#
# Examples:
#   ./notarize/sign_notarize.sh \
#     distribution/macos-apple-silicon/CADENCE_AppleSilicon.app \
#     distribution/macos-apple-silicon/CADENCE-AppleSilicon.dmg  CADENCE
#
#   ./notarize/sign_notarize.sh \
#     distribution/macos-intel/CADENCE_Intel.app \
#     distribution/macos-intel/CADENCE-Intel.dmg  CADENCE

set -euo pipefail

# ============ CONFIG (same for every build) ============
DEVELOPER_ID="Developer ID Application: Kedar Aras (B579TS737G)"   # security find-identity -v -p codesigning
KEYCHAIN_PROFILE="CADENCE_NOTARY"                                  # notarytool store-credentials profile name
ENTITLEMENTS="notarize/entitlements.plist"
# =======================================================

APP="${1:?Usage: $0 <path-to-.app> <output.dmg> [VolumeName]}"
DMG="${2:?Usage: $0 <path-to-.app> <output.dmg> [VolumeName]}"
VOLNAME="${3:-CADENCE}"

[ -d "$APP" ]          || { echo "ERROR: app not found: $APP"; exit 1; }
[ -f "$ENTITLEMENTS" ] || { echo "ERROR: $ENTITLEMENTS not found — run from the project root."; exit 1; }

arch=$(lipo -archs "$APP/Contents/MacOS/$(basename "${APP%.app}")" 2>/dev/null \
       || file "$APP/Contents/MacOS/"* 2>/dev/null | grep -oE 'arm64|x86_64' | head -1)
echo ">> Target: $APP  (arch: ${arch:-unknown})"

echo ">> [1/5] Signing with hardened runtime + entitlements..."
codesign --force --deep --timestamp --options runtime \
         --entitlements "$ENTITLEMENTS" \
         --sign "$DEVELOPER_ID" "$APP"

echo ">> [2/5] Verifying the signature..."
codesign --verify --deep --strict --verbose=2 "$APP"

echo ">> [3/5] Building + signing the DMG..."
STAGE="$(mktemp -d)"
cp -R "$APP" "$STAGE/"
rm -f "$DMG"
hdiutil create -volname "$VOLNAME" -srcfolder "$STAGE" -ov -format UDZO "$DMG"
rm -rf "$STAGE"
codesign --force --timestamp --sign "$DEVELOPER_ID" "$DMG"

echo ">> [4/5] Submitting to Apple for notarization (waits for the result)..."
xcrun notarytool submit "$DMG" --keychain-profile "$KEYCHAIN_PROFILE" --wait

echo ">> [5/5] Stapling the ticket..."
xcrun stapler staple "$DMG"
xcrun stapler staple "$APP" || true

echo
echo ">> Done. Gatekeeper check:"
spctl -a -t open --context context:primary-signature -vv "$DMG" || true
echo ">> Distribute: $DMG"
