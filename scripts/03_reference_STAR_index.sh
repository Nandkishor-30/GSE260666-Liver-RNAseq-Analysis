#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

STAR_BIN="${STAR_BIN:-STAR}"

FASTA="$PROJECT/reference/GRCh38.primary_assembly.genome.fa"
GTF="$PROJECT/reference/gencode.v48.primary_assembly.annotation.gtf"
INDEX_DIR="$PROJECT/star_index"
LOG="$PROJECT/logs/03_STAR_index.log"

mkdir -p "$INDEX_DIR"

"$STAR_BIN" \
  --runMode genomeGenerate \
  --runThreadN "$THREADS" \
  --genomeDir "$INDEX_DIR" \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 150 \
  2>&1 | tee "$LOG"
