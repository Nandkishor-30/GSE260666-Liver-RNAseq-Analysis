dir.create("figures/diagnostics", recursive=TRUE, showWarnings=FALSE)
dir.create("results/diagnostics", recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages({library(DESeq2); library(ggplot2)})
vsd <- readRDS("results/deseq2/vsd_object.rds")
ann <- read.delim("metadata/gencode_v48_gene_annotation.tsv", header=FALSE,
                  col.names=c("gene_id","symbol","biotype"), stringsAsFactors=FALSE)
mk <- c("COL1A1","COL1A2","COL3A1","THBS2","LUM","TIMP1",
        "SPP1","KRT19","EPCAM","CCL2","AKR1B10","CD68")
ids <- ann$gene_id[match(mk, ann$symbol)]; names(ids) <- mk
ids <- ids[!is.na(ids) & ids %in% rownames(vsd)]
cat("markers found:", length(ids), "of", length(mk), "\n")
L <- do.call(rbind, lapply(names(ids), function(g)
  data.frame(gene=g, expr=assay(vsd)[ids[g],], condition=colData(vsd)$condition)))
L$condition <- factor(L$condition, levels=c("Control","NAFL","NASH"))
set.seed(123)
p <- ggplot(L, aes(condition, expr, colour=condition)) +
  geom_boxplot(outlier.shape=NA, alpha=0.3) +
  geom_jitter(width=0.15, size=1.6, alpha=0.85) +
  scale_colour_manual(values = c(Control="#0072B2", NAFL="#E69F00", NASH="#CC79A7")) +
  facet_wrap(~ gene, scales="free_y") +
  labs(y="VST expression", x=NULL, title="Established NAFLD/NASH markers") +
  theme_bw(base_size=9) + theme(legend.position="none")
ggsave("figures/diagnostics/05_marker_panel.png", p,
       width=9, height=7, dpi=300)
write.csv(L, "results/diagnostics/marker_expression.csv", row.names=FALSE)
agg <- aggregate(expr ~ gene + condition, data=L, FUN=mean)
print(reshape(agg, idvar="gene", timevar="condition", direction="wide"))
