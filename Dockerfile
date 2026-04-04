# ── Stage 1: build ──────────────────────────────────────────────────────────
FROM rust:1.87-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    pkg-config libssl-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY Cargo.toml Cargo.lock ./
# Pre-fetch deps (cache layer)
RUN mkdir src && echo 'fn main(){}' > src/main.rs && \
    cargo build --release && rm -rf src

COPY src ./src
RUN touch src/main.rs && cargo build --release

# ── Stage 2: runtime ────────────────────────────────────────────────────────
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/biofastq-a /usr/local/bin/biofastq-a

WORKDIR /data
ENTRYPOINT ["biofastq-a"]
CMD ["--help"]

# ── Usage ────────────────────────────────────────────────────────────────────
# Build:
#   docker build -t biofastq-a .
#
# Run (mount current directory as /data):
#   docker run --rm -v "$PWD":/data biofastq-a sample.fastq --headless
#   docker run --rm -v "$PWD":/data biofastq-a *.fastq.gz --trim --output-dir /data/reports
