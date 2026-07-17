#!/usr/bin/env bash
#
# Runs inside the Linux container (see coverage-local.sh). Mirrors the steps in
# .github/workflows/coverage.yml.
#
set -euo pipefail

# cargo-tarpaulin is installed into this volume-backed root and cached across runs.
TARPAULIN_ROOT=/opt/tarpaulin
export PATH="${TARPAULIN_ROOT}/bin:${PATH}"
export PYO3_PYTHON=python3

echo "==> Installing system dependencies"
apt-get update -qq
apt-get install -y -qq --no-install-recommends \
  pkg-config libssl-dev build-essential binutils python3 python3-dev \
  clang cmake git curl ca-certificates >/dev/null

export RUSTFLAGS="-C link-arg=-fuse-ld=bfd"

if ! command -v cargo-tarpaulin >/dev/null 2>&1; then
  echo "==> Installing cargo-tarpaulin ${TARPAULIN_VERSION}"
  ARCH="$(uname -m)"
  if [[ "${ARCH}" == "x86_64" ]]; then
    mkdir -p "${TARPAULIN_ROOT}/bin"
    curl -sL "https://github.com/xd009642/tarpaulin/releases/download/${TARPAULIN_VERSION}/cargo-tarpaulin-x86_64-unknown-linux-musl.tar.gz" \
      | tar xz -C "${TARPAULIN_ROOT}/bin"
  else
    # No prebuilt musl binary for this arch (e.g. Apple Silicon -> arm64); build it.
    cargo install cargo-tarpaulin --version "${TARPAULIN_VERSION}" --locked --root "${TARPAULIN_ROOT}"
  fi
fi

echo "==> Building binary for CLI tests"
CARGO_BUILD_JOBS=1 cargo build --bins --no-default-features
export CARGO_BIN_EXE_kmerseek="${CARGO_TARGET_DIR}/debug/kmerseek"

echo "==> Generating coverage report"
# --no-dead-code: pyo3 (via sourmash_plugin_branchwater) is pulled in with the
# extension-module feature, so its build script omits the libpython link flags. A
# normal build links fine because --gc-sections drops the unused pyo3 FFI code, but
# tarpaulin defaults to -C link-dead-code which keeps it, leaving PyGen_Type et al.
# unresolved. --no-dead-code turns that off so coverage links like a normal build.
CARGO_BUILD_JOBS=1 cargo tarpaulin --no-dead-code \
  --out Xml --output-dir coverage/ \
  --timeout 120 --skip-clean --no-default-features

echo "==> Done. Report at coverage/cobertura.xml"
