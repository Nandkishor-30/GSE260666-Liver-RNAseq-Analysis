#!/usr/bin/env Rscript

# Read GENCODE annotation
anno <- read.delim(
    "metadata/gencode_v48_gene_annotation.tsv",
    header = FALSE,
    col.names = c("gene_id", "gene_symbol", "gene_type"),
    stringsAsFactors = FALSE
)

comparisons <- c(
    "NAFL_vs_Control",
    "NASH_vs_Control",
    "NASH_vs_NAFL"
)

for (comp in comparisons) {

    # All genes
    infile_all <- paste0(
        "results/deseq2/",
        comp,
        "_all_genes.csv"
    )

    d_all <- read.csv(
        infile_all,
        stringsAsFactors = FALSE
    )

    d_all <- merge(
        d_all,
        anno,
        by = "gene_id",
        all.x = TRUE,
        sort = FALSE
    )

    # Reorder useful columns
    wanted <- c(
        "gene_id",
        "gene_symbol",
        "gene_type",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "significant"
    )

    d_all <- d_all[, wanted]

    write.csv(
        d_all,
        paste0(
            "results/deseq2/",
            comp,
            "_all_genes_annotated.csv"
        ),
        row.names = FALSE
    )

    # Significant DEGs
    d_sig <- subset(
        d_all,
        significant == "Yes"
    )

    d_sig <- d_sig[
        order(d_sig$padj),
    ]

    write.csv(
        d_sig,
        paste0(
            "results/deseq2/",
            comp,
            "_significant_DEGs_annotated.csv"
        ),
        row.names = FALSE
    )

    cat(
        comp,
        ":",
        nrow(d_sig),
        "annotated significant DEGs\n"
    )
}

cat("\nAnnotation completed successfully.\n")
