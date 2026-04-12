#!/bin/bash
set -euxo pipefail

# Disable Rust's own strip step so conda's post-build strip can run cleanly.
export CARGO_PROFILE_RELEASE_STRIP=false

# Use a writable cargo home inside the source tree to avoid $HOME permission
# issues in some conda-build environments.
export CARGO_HOME="${SRC_DIR}/.cargo-home"

cargo install --path . --root "${PREFIX}"
