#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(pheatmap)
    library(DESeq2)
})

dir.create("figures/deseq2", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Volcano plot function
# ------------------------------------------------------------

make_volcano <- function(comp, title_text) {

    f <- paste0(
        "results/deseq2/",
        comp,
        "_all_genes_annotated.csv"
    )

    d <- read.csv(f, stringsAsFactors = FALSE)

    d$category <- "Not significant"

    d$category[
        !is.na(d$padj) &
        d$padj < 0.05 &
        d$log2FoldChange >= 1
    ] <- "Upregulated"

    d$category[
        !is.na(d$padj) &
        d$padj < 0.05 &
        d$log2FoldChange <= -1
    ] <- "Downregulated"

    d$minus_log10_padj <- -log10(d$padj)

    # Avoid Inf values
    finite_vals <- d$minus_log10_padj[is.finite(d$minus_log10_padj)]

    if (length(finite_vals) > 0) {
        max_finite <- max(finite_vals, na.rm = TRUE)
        d$minus_log10_padj[
            is.infinite(d$minus_log10_padj)
        ] <- max_finite + 1
    }

    # Label top significant genes
    sig <- d[
        d$category != "Not significant" &
        !is.na(d$padj),
    ]

    sig <- sig[order(sig$padj), ]

    label_genes <- head(sig, 10)

    p <- ggplot(
        d,
        aes(
            x = log2FoldChange,
            y = minus_log10_padj,
            color = category
        )
    ) +
        geom_point(
            alpha = 0.65,
            size = 1.5
        ) +
        geom_vline(
            xintercept = c(-1, 1),
            linetype = "dashed"
        ) +
        geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed"
        ) +
        geom_text(
            data = label_genes,
            aes(
                label = gene_symbol
            ),
            size = 3,
            vjust = -0.6,
            check_overlap = TRUE,
            show.legend = FALSE
        ) +
        labs(
            title = title_text,
            x = "log2 fold change",
            y = "-log10 adjusted p-value",
            color = "Expression"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(
                hjust = 0.5,
                face = "bold"
            )
        )

    ggsave(
        paste0(
            "figures/deseq2/",
            comp,
            "_volcano.png"
        ),
        p,
        width = 8,
        height = 6,
        dpi = 300
    )
}

# ------------------------------------------------------------
# Three volcano plots
# ------------------------------------------------------------

make_volcano(
    "NAFL_vs_Control",
    "NAFL vs Control"
)

make_volcano(
    "NASH_vs_Control",
    "NASH vs Control"
)

make_volcano(
    "NASH_vs_NAFL",
    "NASH vs NAFL"
)

# ------------------------------------------------------------
# NASH vs Control top-DEG heatmap
# ------------------------------------------------------------

dds <- readRDS(
    "results/deseq2/dds_object.rds"
)

vsd <- readRDS(
    "results/deseq2/vsd_object.rds"
)

res <- read.csv(
    "results/deseq2/NASH_vs_Control_significant_DEGs_annotated.csv",
    stringsAsFactors = FALSE
)

res <- res[
    order(res$padj),
]

# Use top 40 DEGs
top <- head(res, 40)

mat <- assay(vsd)[
    top$gene_id,
    ,
    drop = FALSE
]

# Z-score each gene across samples
mat_z <- t(
    scale(
        t(mat)
    )
)

# Gene labels
gene_labels <- top$gene_symbol

# If gene symbol is missing, retain Ensembl ID
gene_labels[
    is.na(gene_labels) |
    gene_labels == ""
] <- top$gene_id[
    is.na(gene_labels) |
    gene_labels == ""
]

rownames(mat_z) <- make.unique(gene_labels)

metadata <- as.data.frame(
    colData(dds)
)

annotation_col <- data.frame(
    Condition = metadata$condition
)

rownames(annotation_col) <- rownames(metadata)

png(
    "figures/deseq2/NASH_vs_Control_top40_heatmap.png",
    width = 2200,
    height = 2800,
    res = 300
)

pheatmap(
    mat_z,
    annotation_col = annotation_col,
    annotation_colors = list(Condition = c(Control="#0072B2", NAFL="#E69F00", NASH="#CC79A7")),
    show_colnames = TRUE,
    fontsize_row = 7,
    fontsize_col = 7,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Top 40 DEGs: NASH vs Control"
)

dev.off()

cat("Volcano plots and DEG heatmap completed successfully.\n")
