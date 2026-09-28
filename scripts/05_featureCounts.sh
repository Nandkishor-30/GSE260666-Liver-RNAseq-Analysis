#!/bin/bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

FEATURECOUNTS="${FEATURECOUNTS_BIN:-featureCounts}"

mkdir -p counts logs

require_tools "$FEATURECOUNTS" gzip

BAMS=()
declare -A SEEN=()

while IFS= read -r RUN || [[ -n "$RUN" ]]; do
    RUN="${RUN%$'\r'}"
    [[ -z "$RUN" ]] && continue

    if [[ ! "$RUN" =~ ^SRR[0-9]+$ ]]; then
        echo "Invalid run accession: $RUN" >&2
        exit 1
    fi

    if [[ -n "${SEEN[$RUN]:-}" ]]; then
        echo "Duplicate run accession: $RUN" >&2
        exit 1
    fi
    SEEN["$RUN"]=1

    BAM="aligned/${RUN}/${RUN}_Aligned.sortedByCoord.out.bam"
    MARKER="aligned/${RUN}/.alignment_complete"

    if [[ ! -s "$BAM" || ! -f "$MARKER" ]]; then
        echo "Missing BAM or alignment completion marker for $RUN" >&2
        exit 1
    fi

    BAMS+=("$BAM")
done < metadata/SRR_accessions.txt

if [[ "${#BAMS[@]}" -eq 0 ]]; then
    echo "No runs found in metadata/SRR_accessions.txt" >&2
    exit 1
fi

echo "Counting ${#BAMS[@]} BAM files from the accession list."


"$FEATURECOUNTS" \
  -T "$THREADS" \
  -p \
  --countReadPairs \
  -s 2 \
  -t exon \
  -g gene_id \
  -a reference/gencode.v48.primary_assembly.annotation.gtf \
  -o counts/gene_counts.txt \
  "${BAMS[@]}" \
  2>&1 | tee logs/05_featureCounts.log

# compress the count matrix before storing/uploading it - the .summary
# file is tiny and left uncompressed. R reads the .gz directly.
gzip -f counts/gene_counts.txt
