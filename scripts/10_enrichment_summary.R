#!/usr/bin/env Rscript

dir.create("results/final_tables", recursive = TRUE, showWarnings = FALSE)

comparisons <- c("NAFL_vs_Control", "NASH_vs_Control", "NASH_vs_NAFL")
keep_cols <- c("ID", "Description", "NES", "p.adjust", "setSize")

all_summary <- data.frame(
    Comparison = character(),
    Direction = character(),
    ID = character(),
    Description = character(),
    NES = numeric(),
    p.adjust = numeric(),
    setSize = numeric()
)

for (comp in comparisons) {
    infile <- paste0("results/enrichment/", comp, "_GO_BP_GSEA.csv")
    d <- read.csv(infile, stringsAsFactors = FALSE)

    if (nrow(d) == 0L) {
        message(comp, ": no enrichment terms; skipping.")
        next
    }

    if (!all(keep_cols %in% names(d)))
        stop("Missing required columns in: ", infile)

    d <- d[!is.na(d$p.adjust) & d$p.adjust < 0.05 &
           is.finite(d$NES), , drop = FALSE]

    pos <- d[d$NES > 0, , drop = FALSE]
    neg <- d[d$NES < 0, , drop = FALSE]

    pos <- head(pos[order(pos$p.adjust, -pos$NES), , drop = FALSE], 10)
    neg <- head(neg[order(neg$p.adjust, neg$NES), , drop = FALSE], 10)

    tmp <- rbind(pos[, keep_cols, drop = FALSE],
                 neg[, keep_cols, drop = FALSE])

    tmp$Comparison <- rep(comp, nrow(tmp))
    tmp$Direction <- c(
        rep(paste0("Higher in ", sub("_vs_.*", "", comp)), nrow(pos)),
        rep(paste0("Higher in ", sub(".*_vs_", "", comp)), nrow(neg))
    )

    all_summary <- rbind(all_summary, tmp[, names(all_summary), drop = FALSE])
}

write.csv(
    all_summary,
    "results/final_tables/GO_BP_GSEA_top_terms_summary.csv",
    row.names = FALSE
)

infile <- "results/enrichment/NASH_vs_Control_upregulated_GO_BP_ORA.csv"
ora <- read.csv(infile, stringsAsFactors = FALSE)
ora_cols <- c("ID", "Description", "GeneRatio", "BgRatio",
              "FoldEnrichment", "p.adjust", "Count")

if (nrow(ora) == 0L) {
    ora_summary <- data.frame(
        ID = character(), Description = character(),
        GeneRatio = character(), BgRatio = character(),
        FoldEnrichment = numeric(), p.adjust = numeric(), Count = integer()
    )
} else {
    if (!all(ora_cols %in% names(ora)))
        stop("Missing required columns in: ", infile)
    ora_summary <- ora[order(ora$p.adjust), ora_cols, drop = FALSE]
}

write.csv(
    ora_summary,
    "results/final_tables/NASH_vs_Control_upregulated_ORA_summary.csv",
    row.names = FALSE
)

cat("Final enrichment summary tables created successfully.\n")
