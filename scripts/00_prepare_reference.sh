#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
require_tools curl gzip python3
mkdir -p reference metadata
BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48
for FILE in GRCh38.primary_assembly.genome.fa gencode.v48.primary_assembly.annotation.gtf; do
  if [[ ! -s "reference/$FILE" ]]; then
    curl --fail --location --retry 3 "$BASE/$FILE.gz" -o "reference/$FILE.gz.part"
    gzip -t "reference/$FILE.gz.part"
    gzip -dc "reference/$FILE.gz.part" > "reference/$FILE.part"
    mv "reference/$FILE.part" "reference/$FILE"
    rm "reference/$FILE.gz.part"
  fi
done
python3 scripts/extract_gene_annotation.py
