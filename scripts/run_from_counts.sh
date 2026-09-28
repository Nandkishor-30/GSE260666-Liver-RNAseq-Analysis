#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
require_tools Rscript
Rscript -e 'p <- c("DESeq2","ggplot2","pheatmap","clusterProfiler","org.Hs.eg.db","enrichplot","apeglm","ashr"); missing <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))'
for SCRIPT in scripts/0[6-9]_*.R scripts/1[0-4]_*.R; do
  Rscript "$SCRIPT" 2>&1 | tee "logs/$(basename "$SCRIPT" .R).log"
done
