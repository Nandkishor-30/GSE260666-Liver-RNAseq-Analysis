#!/bin/bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"


require_tools fastqc multiqc
mkdir -p qc/raw logs

echo "Starting FastQC: $(date)"

fastqc \
  -t "$THREADS" \
  -o qc/raw \
  fastq/*.fastq.gz \
  2>&1 | tee logs/02_fastqc_all.log

echo "FastQC finished: $(date)"

echo "Starting MultiQC: $(date)"

multiqc qc/raw \
  -o qc \
  --force \
  2>&1 | tee logs/02_multiqc_all.log

echo "MultiQC finished: $(date)"
