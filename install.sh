#!/usr/bin/env bash
# install.sh — build BioFastq-A from source and install the binary
#
# Usage:
#   bash install.sh              # installs to /usr/local/bin (may need sudo)
#   bash install.sh ~/bin        # installs to ~/bin (no sudo)
#   bash install.sh /usr/bin     # custom prefix

set -euo pipefail

BINARY="biofastq-a"
INSTALL_DIR="${1:-/usr/local/bin}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Checks ────────────────────────────────────────────────────────────────────

if ! command -v cargo &>/dev/null; then
    echo "Error: Rust/Cargo not found."
    echo "Install it from https://rustup.rs and re-run this script."
    exit 1
fi

RUST_VERSION=$(rustc --version | awk '{print $2}')
echo "Rust $RUST_VERSION detected."

# ── Build ─────────────────────────────────────────────────────────────────────

echo ""
echo "Building BioFastq-A (release)..."
cd "$REPO_ROOT"
cargo build --release

# ── Install ───────────────────────────────────────────────────────────────────

mkdir -p "$INSTALL_DIR"
DEST="$INSTALL_DIR/$BINARY"

if [[ -w "$INSTALL_DIR" ]]; then
    cp "target/release/$BINARY" "$DEST"
else
    echo "Directory $INSTALL_DIR requires elevated permissions."
    sudo cp "target/release/$BINARY" "$DEST"
fi

chmod 755 "$DEST"

# ── Verify ────────────────────────────────────────────────────────────────────

echo ""
echo "Installed: $DEST"
"$DEST" --version

echo ""
echo "Done!  Try:"
echo "  $BINARY --help"
echo "  $BINARY sample.fastq"
