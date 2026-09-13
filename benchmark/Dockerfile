# mycobench runtime image: Myconductor + mycobench, and nothing else.
#
# Why this image is small and boring
# ----------------------------------
# The bioinformatics tools are NOT baked in here. TB-Profiler and sra-tools
# each ship their own maintained, pinned images, and rebuilding them inside
# this one would mean re-solving their dependency trees badly. So the benchmark
# runs three containers by digest (see docker-compose.yml):
#
#   ncbi/sra-tools      download stage
#   staphb/tbprofiler   profile stage  (carries its own tbdb catalogue)
#   mycobench (this)    interpret / validate / report stages
#
# Myconductor's core has zero runtime dependencies, so this image needs no
# compiler, no conda, and no database download. That is the whole point: the
# interpretation layer stays reproducible and offline-capable even when the
# engines around it are heavy.
#
# Pinned by digest, not by tag: a tag can be repointed.
FROM python:3.12-slim@sha256:5f9c8b4b2d4f2a5cbd0b2e1b0a4e2e4d5e2c8f6e0b1e6c8f5d4a3b2c1d0e9f8a

ARG VCS_REF=unknown
ARG BUILD_DATE=unknown
LABEL org.opencontainers.image.title="mycobench" \
      org.opencontainers.image.description="Real-data validation harness for Myconductor" \
      org.opencontainers.image.source="https://github.com/ubeffiong/myconductor" \
      org.opencontainers.image.revision="$VCS_REF" \
      org.opencontainers.image.created="$BUILD_DATE"

# gzip is used by the download stage to compress fasterq-dump output; it is a
# system tool rather than a Python dependency, and installing it last keeps the
# layer cache intact for the source copy below.
RUN apt-get update \
 && apt-get install -y --no-install-recommends gzip ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/myconductor
COPY . /opt/myconductor
RUN python -m pip install --no-cache-dir --no-deps .

# Non-root: nothing here needs privilege, and the mounted work directories are
# owned by the invoking user.
RUN useradd --create-home --uid 1000 mycobench
USER mycobench

ENV DATA_DIR=/work/data \
    RESULTS_DIR=/work/results \
    LOG_DIR=/work/logs \
    PYTHONUNBUFFERED=1

WORKDIR /work
CMD ["mycobench", "--help"]
