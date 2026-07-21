#!/usr/bin/env bash
#
# Reproduce the GitHub Actions "Coverage" workflow locally in a Linux container.
#
# macOS cannot reproduce the coverage build directly: its linker does dynamic
# symbol lookup by default (so the libpython symbols pyo3 leaves undefined never
# error), and cargo-tarpaulin's coverage engine is Linux-only. This runs the same
# steps as .github/workflows/coverage.yml inside Ubuntu.
#
# Caching: cargo registry + a Linux-only target dir live in named docker volumes,
# so only the first run pays the full compile cost. The host ./target (macOS
# artifacts) is never touched.
#
# Usage:
#   ./scripts/coverage-local.sh              # run coverage, write coverage/cobertura.xml
#   ./scripts/coverage-local.sh --shell      # drop into a shell in the container
#
set -euo pipefail

TARPAULIN_VERSION="0.31.2"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="rust:slim"

if [[ "${1:-}" == "--shell" ]]; then
  ENTRY=(bash)
else
  ENTRY=(bash -euo pipefail /work/scripts/_coverage-inner.sh)
fi

docker run --rm -it \
  -v "${REPO_ROOT}:/work" \
  -v kmerseek-cargo-registry:/usr/local/cargo/registry \
  -v kmerseek-tarpaulin-bin:/opt/tarpaulin \
  -v kmerseek-linux-target:/work/target-linux \
  -e CARGO_TARGET_DIR=/work/target-linux \
  -e CARGO_TERM_COLOR=always \
  -e TARPAULIN_VERSION="${TARPAULIN_VERSION}" \
  -w /work \
  "${IMAGE}" "${ENTRY[@]}"

echo
echo "Coverage report (if generated): ${REPO_ROOT}/coverage/cobertura.xml"
