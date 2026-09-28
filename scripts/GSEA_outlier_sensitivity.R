#!/usr/bin/env Rscript
# Run from the repository root after scripts 06 and 07.
# Original enrichment outputs are preserved.

code <- paste(
  readLines("scripts/09_functional_enrichment.R"),
  collapse = "\n"
)

old <- "!is.na(d$stat) &"
stopifnot(grepl(old, code, fixed = TRUE))

code <- gsub(
  old, "!is.na(d$pvalue) & !is.na(d$stat) &",
  code, fixed = TRUE
)
code <- gsub(
  "results/enrichment",
  "results/enrichment_outlier_sensitivity",
  code, fixed = TRUE
)
code <- gsub(
  "figures/enrichment",
  "figures/enrichment_outlier_sensitivity",
  code, fixed = TRUE
)

eval(parse(text = code))

comparisons <- c(
  "NAFL_vs_Control", "NASH_vs_Control", "NASH_vs_NAFL"
)

summary <- do.call(rbind, lapply(comparisons, function(comp) {
  filename <- paste0(comp, "_GO_BP_GSEA.csv")
  old <- read.csv(file.path("results/enrichment", filename))
  new <- read.csv(file.path(
    "results/enrichment_outlier_sensitivity", filename
  ))
  shared <- merge(old, new, by = "ID", suffixes = c("_old", "_new"))

  data.frame(
    comparison = comp,
    original_terms = nrow(old),
    sensitivity_terms = nrow(new),
    shared_terms = nrow(shared),
    retained_percent = round(100 * nrow(shared) / nrow(old), 1),
    lost_terms = sum(!old$ID %in% new$ID),
    gained_terms = sum(!new$ID %in% old$ID),
    direction_reversals = sum(
      sign(shared$NES_old) != sign(shared$NES_new)
    )
  )
}))

dir.create("validation", showWarnings = FALSE)
write.csv(summary, "validation/GSEA_sensitivity_summary.csv",
          row.names = FALSE)
capture.output(sessionInfo(),
               file = "validation/GSEA_sensitivity_sessionInfo.txt")
print(summary, row.names = FALSE)
