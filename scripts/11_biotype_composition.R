dir.create("results/diagnostics", recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(DESeq2))
dds <- readRDS("results/deseq2/dds_object.rds")
cat("dds:", nrow(dds), "genes x", ncol(dds), "samples\n")
print(table(colData(dds)$condition))
ann <- read.delim("metadata/gencode_v48_gene_annotation.tsv", header=FALSE,
                  col.names=c("gene_id","symbol","biotype"), stringsAsFactors=FALSE)
raw <- counts(dds, normalized=FALSE)
bt <- ann$biotype[match(rownames(raw), ann$gene_id)]
bt[is.na(bt)] <- "unknown"
agg <- rowsum(raw, group=bt)
pct <- round(100*t(t(agg)/colSums(agg)), 3)
keep <- order(rowMeans(pct), decreasing=TRUE)[1:12]
cat("\nTOP BIOTYPES (mean % of counts):\n")
print(round(rowMeans(pct)[keep], 2))
rr <- grep("rRNA", rownames(pct), value=TRUE)
cat("\nrRNA share per sample (%):\n")
print(round(colSums(pct[rr, , drop=FALSE]), 3))
write.csv(as.data.frame(pct[keep, ]), "results/diagnostics/biotype_percent.csv")
