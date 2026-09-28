#!/bin/bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"


require_tools prefetch vdb-validate fasterq-dump pigz gzip
mkdir -p sra fastq tmp logs

while read -r RUN
do
    echo "=============================================="
    echo "Processing $RUN"
    echo "Started: $(date)"
    echo "=============================================="

    # Skip runs for which both compressed FASTQs already exist
    if [[ -f "fastq/${RUN}_1.fastq.gz" && \
          -f "fastq/${RUN}_2.fastq.gz" ]]; then
        gzip -t "fastq/${RUN}_1.fastq.gz" "fastq/${RUN}_2.fastq.gz"
        echo "$RUN already completed. Skipping."
        echo
        continue
    fi

    # Download only if the SRA archive is absent
    if [[ ! -f "sra/${RUN}/${RUN}.sra" ]]; then
        echo "Downloading $RUN..."

        prefetch "$RUN" \
            --output-directory sra \
            2>&1 | tee "logs/${RUN}_prefetch.log"
    else
        echo "Existing SRA archive found for $RUN."
    fi

    echo "Validating $RUN..."

    vdb-validate "sra/${RUN}/${RUN}.sra" \
        2>&1 | tee "logs/${RUN}_validate.log"

    # Convert in a fresh staging directory so interrupted extraction is never reused.
    if [[ ! -s "fastq/${RUN}_1.fastq.gz" || ! -s "fastq/${RUN}_2.fastq.gz" ]]; then
        STAGE=$(mktemp -d "tmp/${RUN}.XXXXXX")
        fasterq-dump "sra/${RUN}/${RUN}.sra" --split-files \
          --threads "$THREADS" --temp "$STAGE/work" --outdir "$STAGE" \
          2>&1 | tee "logs/${RUN}_fasterq.log"
        test -s "$STAGE/${RUN}_1.fastq"
        test -s "$STAGE/${RUN}_2.fastq"
        pigz -p "$THREADS" "$STAGE/${RUN}_1.fastq" "$STAGE/${RUN}_2.fastq"
        gzip -t "$STAGE/${RUN}_1.fastq.gz" "$STAGE/${RUN}_2.fastq.gz"
        mv -f "$STAGE/${RUN}_1.fastq.gz" "$STAGE/${RUN}_2.fastq.gz" fastq/
        rm -rf "$STAGE"
    fi

    # Final integrity check
    if [[ -f "fastq/${RUN}_1.fastq.gz" && \
          -f "fastq/${RUN}_2.fastq.gz" ]]; then

        echo "Finished $RUN successfully"
        ls -lh \
            "fastq/${RUN}_1.fastq.gz" \
            "fastq/${RUN}_2.fastq.gz"

    else
        echo "ERROR: final FASTQ pair missing for $RUN"
        exit 1
    fi

    echo "Completed: $(date)"
    echo

done < metadata/SRR_accessions.txt
