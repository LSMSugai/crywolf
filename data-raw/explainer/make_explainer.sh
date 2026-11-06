#!/usr/bin/env bash
# Draw the panels of vignettes/articles/how-crywolf-works.Rmd into vignettes/articles/figures/.
# Needs: python3, node, and Playwright with Chromium (npm install playwright && npx playwright install chromium).
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -d fontpkg/node_modules/@fontsource ]; then
  mkdir -p fontpkg && (cd fontpkg && npm init -y >/dev/null && npm install @fontsource/figtree@5 @fontsource/ibm-plex-mono@5 >/dev/null)
fi
mkdir -p ../../vignettes/articles/figures
python3 build.py
node render.js
echo "panels written to vignettes/articles/figures/how-*.png"
