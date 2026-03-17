# ── builder ──────────────────────────────────────────────────────────────────
FROM rust:slim AS builder

# RocksDB needs clang/llvm; git is required for the sourmash_plugin_branchwater
# git dependency in Cargo.toml
RUN apt-get update && apt-get install -y --no-install-recommends \
    clang \
    cmake \
    git \
    libclang-dev \
    python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

COPY Cargo.toml Cargo.lock rust-toolchain.toml ./
COPY src/ src/
COPY benches/ benches/
COPY examples/ examples/

RUN cargo build --release --no-default-features --bin kmerseek

# ── runtime ──────────────────────────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/kmerseek /usr/local/bin/kmerseek

WORKDIR /data

ENTRYPOINT ["kmerseek"]
