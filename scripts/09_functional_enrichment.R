#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
})

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/enrichment", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper: remove GENCODE version suffix from Ensembl IDs
# ------------------------------------------------------------

strip_version <- function(x) {
    sub("\\.[0-9]+$", "", x)
}

# ------------------------------------------------------------
# Helper: convert Ensembl IDs to Entrez IDs
# ------------------------------------------------------------

map_to_entrez <- function(ensembl_ids) {

    clean_ids <- unique(strip_version(ensembl_ids))

    mapped <- bitr(
        clean_ids,
        fromType = "ENSEMBL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Hs.eg.db
    )

    mapped
}

# ------------------------------------------------------------
# 1. NASH vs Control ORA
# ------------------------------------------------------------

all_nash <- read.csv(
    "results/deseq2/NASH_vs_Control_all_genes_annotated.csv",
    stringsAsFactors = FALSE
)

sig_nash <- read.csv(
    "results/deseq2/NASH_vs_Control_significant_DEGs_annotated.csv",
    stringsAsFactors = FALSE
)

# Background/universe = genes tested by DESeq2 and mappable to Entrez
universe_map <- map_to_entrez(all_nash$gene_id)
universe_entrez <- unique(universe_map$ENTREZID)

# Significant genes
sig_map <- map_to_entrez(sig_nash$gene_id)
sig_entrez <- unique(sig_map$ENTREZID)

# Upregulated
up_ids <- sig_nash$gene_id[
    sig_nash$log2FoldChange >= 1
]

up_map <- map_to_entrez(up_ids)
up_entrez <- unique(up_map$ENTREZID)

# Downregulated
down_ids <- sig_nash$gene_id[
    sig_nash$log2FoldChange <= -1
]

down_map <- map_to_entrez(down_ids)
down_entrez <- unique(down_map$ENTREZID)

cat("NASH vs Control significant DEGs:", nrow(sig_nash), "\n")
cat("Mapped significant DEGs:", length(sig_entrez), "\n")
cat("Mapped upregulated DEGs:", length(up_entrez), "\n")
cat("Mapped downregulated DEGs:", length(down_entrez), "\n")
cat("Mapped background genes:", length(universe_entrez), "\n")

# ------------------------------------------------------------
# GO Biological Process ORA function
# ------------------------------------------------------------

run_enrichGO <- function(gene_ids, label) {

    ego <- enrichGO(
        gene = gene_ids,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
    )

    result_df <- as.data.frame(ego)

    write.csv(
        result_df,
        paste0(
            "results/enrichment/",
            label,
            "_GO_BP_ORA.csv"
        ),
        row.names = FALSE
    )

    cat(
        label,
        "- significant GO BP terms:",
        nrow(result_df),
        "\n"
    )

    if (nrow(result_df) > 0) {

        p <- dotplot(
            ego,
            showCategory = min(15, nrow(result_df))
        ) +
            ggtitle(
                paste0(
                    label,
                    " - GO Biological Process"
                )
            )

        ggsave(
            paste0(
                "figures/enrichment/",
                label,
                "_GO_BP_ORA_dotplot.png"
            ),
            p,
            width = 9,
            height = 7,
            dpi = 300
        )
    }

    return(ego)
}

# ------------------------------------------------------------
# Run ORA
# ------------------------------------------------------------

ego_all <- run_enrichGO(
    sig_entrez,
    "NASH_vs_Control_all_DEGs"
)

ego_up <- run_enrichGO(
    up_entrez,
    "NASH_vs_Control_upregulated"
)

ego_down <- run_enrichGO(
    down_entrez,
    "NASH_vs_Control_downregulated"
)

# ------------------------------------------------------------
# 2. GSEA helper
# ------------------------------------------------------------

make_ranked_gene_list <- function(filename) {

    d <- read.csv(
        filename,
        stringsAsFactors = FALSE
    )

    # Keep genes with a valid DESeq2 Wald statistic
    d <- d[
        !is.na(d$stat) &
        is.finite(d$stat),
    ]

    d$ENSEMBL_clean <- strip_version(d$gene_id)

    mapping <- bitr(
        unique(d$ENSEMBL_clean),
        fromType = "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
    )

    merged <- merge(
        d,
        mapping,
        by.x = "ENSEMBL_clean",
        by.y = "ENSEMBL"
    )

    # If multiple Ensembl IDs map to the same Entrez ID,
    # retain the gene with the largest absolute Wald statistic.
    merged <- merged[
        order(abs(merged$stat), decreasing = TRUE),
    ]

    merged <- merged[
        !duplicated(merged$ENTREZID),
    ]

    gene_list <- merged$stat
    names(gene_list) <- merged$ENTREZID

    gene_list <- sort(
        gene_list,
        decreasing = TRUE
    )

    return(gene_list)
}

# ------------------------------------------------------------
# GO-BP GSEA function
# ------------------------------------------------------------

run_gsea <- function(comp) {

    infile <- paste0(
        "results/deseq2/",
        comp,
        "_all_genes_annotated.csv"
    )

    gene_list <- make_ranked_gene_list(infile)

    cat(
        comp,
        "- genes in ranked GSEA list:",
        length(gene_list),
        "\n"
    )

    set.seed(123)

    gsea <- gseGO(
        geneList = gene_list,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        eps = 0,
        seed = TRUE,
        verbose = FALSE
    )

    gsea_df <- as.data.frame(gsea)

    write.csv(
        gsea_df,
        paste0(
            "results/enrichment/",
            comp,
            "_GO_BP_GSEA.csv"
        ),
        row.names = FALSE
    )

    cat(
        comp,
        "- significant GSEA GO BP terms:",
        nrow(gsea_df),
        "\n"
    )

    if (nrow(gsea_df) > 0) {

        # NOTE: scales = "free_y", space = "free_y" lets each facet panel
        # size its own y-axis instead of sharing one. That alone is not
        # enough when GO term names are long and wrap onto two lines —
        # with ~15 categories per panel there isn't enough vertical room
        # per label, so wrapped text from adjacent terms still collides.
        # Fixing this also requires: (1) a taller image so ggplot has
        # more vertical space to spread labels apart, and (2) smaller
        # axis text so each wrapped label takes less vertical room.
        n_cat <- min(15, nrow(gsea_df))

        p <- dotplot(
            gsea,
            showCategory = n_cat,
            split = ".sign"
        ) +
            facet_grid(. ~ .sign, scales = "free_y", space = "free_y") +
            ggtitle(
                paste0(
                    comp,
                    " - GO BP GSEA"
                )
            ) +
            theme(
                axis.text.y = element_text(size = 8),
                strip.text = element_text(size = 10)
            )

        # Scale image height to the number of categories so labels have
        # room to breathe regardless of how many terms are plotted.
        plot_height <- max(7, n_cat * 0.55)

        ggsave(
            paste0(
                "figures/enrichment/",
                comp,
                "_GO_BP_GSEA_dotplot.png"
            ),
            p,
            width = 12,
            height = plot_height,
            dpi = 300
        )
    }

    return(gsea)
}

# ------------------------------------------------------------
# Run GSEA for all three comparisons
# ------------------------------------------------------------

gsea_nafl_control <- run_gsea(
    "NAFL_vs_Control"
)

gsea_nash_control <- run_gsea(
    "NASH_vs_Control"
)

gsea_nash_nafl <- run_gsea(
    "NASH_vs_NAFL"
)

# ------------------------------------------------------------
# Save enrichment objects for reproducibility
# ------------------------------------------------------------

saveRDS(
    list(
        ORA_all = ego_all,
        ORA_up = ego_up,
        ORA_down = ego_down,
        GSEA_NAFL_vs_Control = gsea_nafl_control,
        GSEA_NASH_vs_Control = gsea_nash_control,
        GSEA_NASH_vs_NAFL = gsea_nash_nafl
    ),
    "results/enrichment/enrichment_objects.rds"
)

cat("\nFunctional enrichment analysis completed successfully.\n")
