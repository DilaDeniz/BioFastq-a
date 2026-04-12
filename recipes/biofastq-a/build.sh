#!/bin/bash
set -euxo pipefail

export CARGO_PROFILE_RELEASE_STRIP=false
cargo build --release
install -m755 target/release/biofastq-a "$PREFIX/bin/biofastq-a"
