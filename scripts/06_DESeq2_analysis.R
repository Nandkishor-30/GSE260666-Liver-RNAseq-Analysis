#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/deseq2", recursive = TRUE, showWarnings = FALSE)

# Read featureCounts output
fc <- read.delim(
    "counts/gene_counts.txt.gz",
    comment.char = "#",
    check.names = FALSE
)

cat("Genes in featureCounts file:", nrow(fc), "\n")
cat("Columns in featureCounts file:", ncol(fc), "\n")

count_matrix <- fc[, 7:ncol(fc)]
rownames(count_matrix) <- fc$Geneid

colnames(count_matrix) <- basename(colnames(count_matrix))
colnames(count_matrix) <- sub(
    "_Aligned.sortedByCoord.out.bam$",
    "",
    colnames(count_matrix)
)

count_matrix <- as.matrix(count_matrix)

# Validate raw counts before integer conversion.
if (!is.numeric(count_matrix))
    stop("Count matrix must contain numeric values")
if (anyNA(count_matrix) || any(!is.finite(count_matrix)))
    stop("Count matrix contains missing or non-finite values")
if (any(count_matrix < 0))
    stop("Count matrix contains negative values")
if (any(count_matrix != floor(count_matrix)))
    stop("Count matrix contains fractional values")
if (any(count_matrix > .Machine$integer.max))
    stop("Count matrix exceeds R's integer range")

storage.mode(count_matrix) <- "integer"

# Read metadata
metadata <- read.delim(
    "metadata/sample_metadata.tsv",
    header = TRUE,
    stringsAsFactors = FALSE
)

if (anyDuplicated(metadata$sample) || anyDuplicated(colnames(count_matrix)))
    stop("Duplicate sample identifiers")
if (!setequal(metadata$sample, colnames(count_matrix)))
    stop("Metadata and count matrix must contain exactly the same samples")
if (anyNA(metadata$condition) || !all(metadata$condition %in% c("Control", "NAFL", "NASH")))
    stop("Missing or invalid condition")
rownames(metadata) <- metadata$sample

metadata$condition <- factor(
    metadata$condition,
    levels = c("Control", "NAFL", "NASH")
)

# Verify sample matching
if (!all(colnames(count_matrix) %in% rownames(metadata))) {
    stop("Some count matrix samples are missing from metadata.")
}

metadata <- metadata[colnames(count_matrix), , drop = FALSE]

if (!identical(colnames(count_matrix), rownames(metadata))) {
    stop("Count matrix and metadata sample order does not match.")
}

cat("\nSample matching successful.\n")
print(table(metadata$condition))

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ condition
)

# Low-expression filter
keep <- rowSums(counts(dds) >= 10) >= 4

cat("\nGenes before filtering:", nrow(dds), "\n")

dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")

# Run DESeq2
dds <- DESeq(dds)

saveRDS(
    dds,
    "results/deseq2/dds_object.rds"
)

# Normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

write.csv(
    normalized_counts,
    "results/deseq2/normalized_counts.csv"
)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

saveRDS(
    vsd,
    "results/deseq2/vsd_object.rds"
)

# PCA
pcaData <- plotPCA(
    vsd,
    intgroup = "condition",
    returnData = TRUE
)

percentVar <- round(
    100 * attr(pcaData, "percentVar")
)

p <- ggplot(
    pcaData,
    aes(
        x = PC1,
        y = PC2,
        color = condition,
        label = name
    )
) +
    geom_point(size = 4) +
    geom_text(
        vjust = -0.8,
        size = 3,
        show.legend = FALSE
    ) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of GSE260666 liver RNA-seq samples") +
    theme_bw()

ggsave(
    "figures/deseq2/PCA_conditions.png",
    p,
    width = 8,
    height = 6,
    dpi = 300
)

write.csv(
    pcaData,
    "results/deseq2/PCA_coordinates.csv"
)

# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- metadata$sample
colnames(sampleDistMatrix) <- metadata$sample

annotation <- data.frame(
    Condition = metadata$condition
)

rownames(annotation) <- metadata$sample

png(
    "figures/deseq2/sample_distance_heatmap.png",
    width = 2200,
    height = 2000,
    res = 300
)

pheatmap(
    sampleDistMatrix,
    annotation_col = annotation,
    annotation_row = annotation,
    main = "Sample-to-sample distance"
)

dev.off()

# Differential expression helper
run_comparison <- function(dds, numerator, denominator, filename) {

    res <- results(
        dds,
        contrast = c(
            "condition",
            numerator,
            denominator
        ),
        alpha = 0.05
    )

    res <- res[order(res$padj), ]

    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)

    res_df$significant <- ifelse(
        !is.na(res_df$padj) &
        res_df$padj < 0.05 &
        abs(res_df$log2FoldChange) >= 1,
        "Yes",
        "No"
    )

    write.csv(
        res_df,
        paste0(
            "results/deseq2/",
            filename,
            "_all_genes.csv"
        ),
        row.names = FALSE
    )

    significant <- subset(
        res_df,
        significant == "Yes"
    )

    write.csv(
        significant,
        paste0(
            "results/deseq2/",
            filename,
            "_significant_DEGs.csv"
        ),
        row.names = FALSE
    )

    cat(
        "\n",
        numerator,
        "vs",
        denominator,
        "\n"
    )

    cat(
        "Significant DEGs:",
        nrow(significant),
        "\n"
    )

    cat(
        "Upregulated:",
        sum(significant$log2FoldChange >= 1),
        "\n"
    )

    cat(
        "Downregulated:",
        sum(significant$log2FoldChange <= -1),
        "\n"
    )
}

run_comparison(
    dds,
    "NAFL",
    "Control",
    "NAFL_vs_Control"
)

run_comparison(
    dds,
    "NASH",
    "Control",
    "NASH_vs_Control"
)

run_comparison(
    dds,
    "NASH",
    "NAFL",
    "NASH_vs_NAFL"
)

sink(
    "results/deseq2/R_sessionInfo.txt"
)

sessionInfo()

sink()

cat("\nDESeq2 analysis completed successfully.\n")
