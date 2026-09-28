#!/bin/bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

STAR="${STAR_BIN:-STAR}"

mkdir -p aligned logs

while read -r RUN
do
    echo "=============================================="
    echo "Aligning $RUN"
    echo "Started: $(date)"
    echo "=============================================="

    OUTDIR="aligned/${RUN}"
    mkdir -p "$OUTDIR"

    if [[ -s "${OUTDIR}/${RUN}_Aligned.sortedByCoord.out.bam" && -f "${OUTDIR}/.alignment_complete" ]]; then
        echo "BAM already exists for $RUN - skipping"
        continue
    fi

    "$STAR" \
        --runThreadN "$THREADS" \
        --genomeDir star_index \
        --readFilesIn \
            "fastq/${RUN}_1.fastq.gz" \
            "fastq/${RUN}_2.fastq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${OUTDIR}/${RUN}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM \
        --quantMode GeneCounts \
        --outSAMunmapped Within \
        --outSAMstrandField intronMotif

    touch "${OUTDIR}/.alignment_complete"
    echo "Finished $RUN"
    echo "Completed: $(date)"

done < metadata/SRR_accessions.txt
