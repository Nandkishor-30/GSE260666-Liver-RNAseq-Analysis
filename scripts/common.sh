#!/usr/bin/env bash
# Shared portable settings; activate dependencies before running the pipeline.
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT"
THREADS="${THREADS:-8}"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 1; }
mkdir -p logs
require_tools() {
  for tool in "$@"; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Missing executable: $tool" >&2; exit 1; }
  done
}
