#!/bin/bash
# Sign, notarize, and staple the CADENCE standalone app and a distributable DMG.
# Prereqs: Developer ID Application cert installed, and a stored notarytool
# credential profile (see README in this folder / the chat instructions).
#
# Usage:  ./notarize/sign_notarize.sh
# Run from the project root (/Users/kedararas/Desktop/Cadence2026).

set -euo pipefail

# ============ FILL THESE IN ============
DEVELOPER_ID="Developer ID Application: Kedar Aras (B579TS737G)"   # exact string from: security find-identity -v -p codesigning
KEYCHAIN_PROFILE="CADENCE_NOTARY"                              # the profile name you stored with notarytool store-credentials
# =======================================

APP="CADENCE_Desktop_App/build/CADENCE.app"
ENTITLEMENTS="notarize/entitlements.plist"
DMG="CADENCE_Desktop_App/CADENCE.dmg"
VOLNAME="CADENCE"

[ -d "$APP" ] || { echo "ERROR: $APP not found (run from project root)."; exit 1; }

echo ">> [1/5] Signing the app with hardened runtime + entitlements..."
# --deep signs nested binaries too; works for most MATLAB standalone apps.
# If notarization later reports an unsigned nested binary, switch to an
# inside-out sign (sign each Mach-O under Contents/ first, then the bundle).
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
codesign -dv --verbose=4 "$APP" 2>&1 | grep -E "Authority|TeamIdentifier|Timestamp" || true
echo ">> Distribute: $DMG"
