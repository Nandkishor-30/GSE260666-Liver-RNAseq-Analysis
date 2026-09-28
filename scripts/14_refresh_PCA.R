suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

if (!requireNamespace("ggrepel", quietly = TRUE))
  stop("ggrepel is missing; install it before continuing.")

group_colors <- c(
  Control = "#0072B2",
  NAFL = "#E69F00",
  NASH = "#CC79A7"
)

vsd <- readRDS("results/deseq2/vsd_object.rds")
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
variance <- round(100 * attr(pca, "percentVar"), 1)

p <- ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3.5) +
  ggrepel::geom_text_repel(
    aes(label = name),
    seed = 123,
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(
    title = "PCA of human liver RNA-seq samples",
    subtitle = "Control (n = 6), NAFL (n = 6), NASH (n = 4)",
    x = paste0("PC1: ", variance[1], "% variance"),
    y = paste0("PC2: ", variance[2], "% variance"),
    color = "Condition"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  "figures/deseq2/PCA_conditions.png",
  p, width = 10, height = 7, dpi = 300
)

cat("PCA figure updated successfully.\n")

# Refresh sample-distance heatmap
sample_dist <- as.matrix(dist(t(assay(vsd))))
annotation <- data.frame(Condition = colData(vsd)$condition)
rownames(annotation) <- colnames(vsd)

png("figures/deseq2/sample_distance_heatmap.png",
    width = 2600, height = 2400, res = 300)
pheatmap::pheatmap(
  sample_dist,
  annotation_col = annotation,
  annotation_row = annotation,
  annotation_colors = list(Condition = group_colors),
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Sample-to-sample distance"
)
dev.off()
cat("Sample-distance heatmap updated successfully.\n")
