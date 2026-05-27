#!/usr/bin/env bash
# Build the BioFastq-A WebAssembly module and copy output to www/pkg/.
# Prerequisites: wasm-pack  (cargo install wasm-pack)
#                           or  curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "==> Building WebAssembly module (release)…"
wasm-pack build --target web --out-dir www/pkg --release

echo ""
echo "==> Build complete. Output: wasm/www/pkg/"
echo ""
echo "To serve locally:"
echo "  cd wasm/www"
echo "  python3 -m http.server 8080"
echo "  # then open http://localhost:8080"
echo ""
echo "To deploy: copy the entire wasm/www/ directory to any static host"
echo "  (GitHub Pages, Netlify, Vercel, S3, etc.)"
