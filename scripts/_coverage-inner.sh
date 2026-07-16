#!/usr/bin/env bash
#
# Runs inside the Linux container (see coverage-local.sh). Mirrors the steps in
# .github/workflows/coverage.yml.
#
set -euo pipefail

echo "==> Installing system dependencies"
apt-get update -qq
apt-get install -y -qq --no-install-recommends \
  pkg-config libssl-dev build-essential binutils python3-dev \
  clang cmake git curl ca-certificates >/dev/null

PYTHON_LIB_DIR=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
PYTHON_LDVERSION=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('LDVERSION'))")
export LD_LIBRARY_PATH="${PYTHON_LIB_DIR}:${LD_LIBRARY_PATH:-}"

# Same reasoning as the workflow: pyo3 (via sourmash_plugin_branchwater) uses the
# extension-module feature, so its build script omits the libpython link flags.
# tarpaulin builds with -C link-dead-code, keeping the unused pyo3 FFI code, so we
# must link libpython explicitly. --no-as-needed is required: rustc emits -lpython
# before the pyo3 objects, and the linker's default --as-needed would drop it.
PYTHON_LINK_FLAGS="-C link-arg=-L${PYTHON_LIB_DIR} -C link-arg=-Wl,--no-as-needed -C link-arg=-lpython${PYTHON_LDVERSION} -C link-arg=-Wl,--as-needed"
export RUSTFLAGS="-C link-arg=-fuse-ld=bfd ${PYTHON_LINK_FLAGS}"

echo "==> RUSTFLAGS=${RUSTFLAGS}"

if ! command -v cargo-tarpaulin >/dev/null 2>&1; then
  echo "==> Installing cargo-tarpaulin ${TARPAULIN_VERSION}"
  ARCH="$(uname -m)"
  if [[ "${ARCH}" == "x86_64" ]]; then
    curl -sL "https://github.com/xd009642/tarpaulin/releases/download/${TARPAULIN_VERSION}/cargo-tarpaulin-x86_64-unknown-linux-musl.tar.gz" \
      | tar xz -C /usr/local/cargo/bin
  else
    # No prebuilt musl binary for this arch (e.g. Apple Silicon -> arm64); build it.
    cargo install cargo-tarpaulin --version "${TARPAULIN_VERSION}" --locked
  fi
fi

echo "==> Building binary for CLI tests"
CARGO_BUILD_JOBS=1 cargo build --bins --no-default-features
export CARGO_BIN_EXE_kmerseek="${CARGO_TARGET_DIR}/debug/kmerseek"

echo "==> Generating coverage report"
CARGO_BUILD_JOBS=1 cargo tarpaulin \
  --out Xml --output-dir coverage/ \
  --timeout 120 --skip-clean --no-default-features

echo "==> Done. Report at coverage/cobertura.xml"
