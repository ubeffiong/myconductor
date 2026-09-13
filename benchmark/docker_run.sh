#!/usr/bin/env bash
# Run one benchmark service, recording the image's immutable identity.
#
# A published number has to be traceable to the exact image that produced it.
# A tag is not that: it can be repointed. So this script resolves the running
# image to its repo digest (falling back to the local image id when the image
# was built locally and never pushed) and passes it into the container, where
# the run manifest records it.
#
# Usage:
#   benchmark/docker_run.sh sra        prefetch -O /work/data/ng_SRR5535815 SRR5535815
#   benchmark/docker_run.sh tbprofiler tb-profiler version
#   benchmark/docker_run.sh mycobench  mycobench run --cohort nigeria-v1
set -euo pipefail

SERVICE="${1:-}"
shift || true

case "$SERVICE" in
  sra)
    IMAGE="ncbi/sra-tools:3.4.1@sha256:1bf9aa259adc315df1cf9d6f8e9d5de7146f1d167e4b7abdc7a61438423da795"
    ;;
  tbprofiler)
    IMAGE="staphb/tbprofiler:6.7.0@sha256:c3785d4267c21b3ea6c083a49260271753054954750f39f8cea38056d7390eee"
    ;;
  mycobench)
    IMAGE="${MYCOBENCH_IMAGE:-mycobench:local}"
    ;;
  *)
    echo "usage: $0 {sra|tbprofiler|mycobench} COMMAND..." >&2
    echo "  sra         download stage (ncbi/sra-tools)" >&2
    echo "  tbprofiler  profile stage (staphb/tbprofiler)" >&2
    echo "  mycobench   interpret / validate / report stages" >&2
    exit 2
    ;;
esac

if ! docker info >/dev/null 2>&1; then
  echo "ERROR: the Docker daemon is not reachable." >&2
  echo "  Start Docker Desktop (or dockerd) and try again. Nothing was run." >&2
  exit 1
fi

# Resolve an immutable identity for the manifest. A digest-pinned reference
# already is one; a locally built image has no repo digest, so fall back to its
# content-addressed image id rather than recording a bare tag.
digest="$(docker image inspect "$IMAGE" --format '{{if .RepoDigests}}{{index .RepoDigests 0}}{{else}}{{.Id}}{{end}}' 2>/dev/null || true)"
if [[ -z "$digest" || "$digest" == "<no value>" ]]; then
  digest="$IMAGE (not yet pulled; digest unresolved)"
fi

HOST_ROOT="${MYCOBENCH_HOST_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
mkdir -p "$HOST_ROOT/data" "$HOST_ROOT/results" "$HOST_ROOT/logs"

exec docker run --rm \
  -v "$HOST_ROOT/cohorts:/work/cohorts:ro" \
  -v "$HOST_ROOT/data:/work/data" \
  -v "$HOST_ROOT/results:/work/results" \
  -v "$HOST_ROOT/logs:/work/logs" \
  -e DATA_DIR=/work/data \
  -e RESULTS_DIR=/work/results \
  -e LOG_DIR=/work/logs \
  -e "NCBI_EMAIL=${NCBI_EMAIL:-}" \
  -e "NCBI_API_KEY=${NCBI_API_KEY:-}" \
  -e "CONTAINER_IMAGE=$IMAGE" \
  -e "CONTAINER_IMAGE_DIGEST=$digest" \
  -w /work \
  "$IMAGE" "$@"
