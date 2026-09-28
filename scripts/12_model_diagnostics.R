dir.create("figures/diagnostics", recursive=TRUE, showWarnings=FALSE)
dir.create("results/diagnostics", recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(DESeq2))
dds <- readRDS("results/deseq2/dds_object.rds")
FIG <- "figures/diagnostics/"
TAB <- "results/diagnostics/"
png(paste0(FIG,"01_dispersion.png"), width=1800, height=1500, res=300)
plotDispEsts(dds, main="DESeq2 dispersion estimates"); dev.off()
cs <- list(NAFL_vs_Control=c("NAFL","Control"), NASH_vs_Control=c("NASH","Control"),
           NASH_vs_NAFL=c("NASH","NAFL")); out <- list()
for (nm in names(cs)) {
  a <- cs[[nm]][1]; b <- cs[[nm]][2]
  r <- results(dds, contrast=c("condition",a,b), alpha=0.05)
  cf <- paste0("condition_",a,"_vs_",b)
  method <- if (cf %in% resultsNames(dds)) "apeglm" else "ashr"
  if (!requireNamespace(method, quietly=TRUE))
    stop("Install the required shrinkage package: ", method)
  # Fail explicitly: unshrunken values must never be labelled as shrunken.
  s <- if (method == "apeglm") lfcShrink(dds, coef=cf, type=method) else
    lfcShrink(dds, contrast=c("condition",a,b), res=r, type=method)
  d <- as.data.frame(r); d$gene_id <- rownames(d)
  d$lfc_shrunk <- s$log2FoldChange[match(rownames(d), rownames(s))]
  write.csv(d[order(d$padj),], paste0(TAB,"shrunken_",nm,".csv"), row.names=FALSE)
  png(paste0(FIG,"02_pval_",nm,".png"), width=1700, height=1300, res=300)
  hist(d$pvalue, breaks=40, col="grey40", border="white", xlab="p-value", main=nm); dev.off()
  png(paste0(FIG,"03_MA_",nm,".png"), width=1800, height=1500, res=300)
  L <- ifelse(is.na(d$lfc_shrunk), d$log2FoldChange, d$lfc_shrunk)
  plot(d$baseMean, L, log="x", pch=16, cex=0.35, ylim=c(-6,6),
       col=ifelse(!is.na(d$padj) & d$padj<0.05, "red", "grey65"),
       xlab="mean normalised count", ylab="shrunken log2FC", main=nm)
  abline(h=0, col="blue"); abline(h=c(-1,1), lty=2, col="blue"); dev.off()
  k <- !is.na(d$padj)
  out[[nm]] <- data.frame(contrast=nm, tested=sum(k), padj05=sum(d$padj<0.05,na.rm=TRUE),
    sig=sum(k & d$padj<0.05 & abs(L)>=1, na.rm=TRUE))
}
res <- do.call(rbind,out); print(res)
write.csv(res, paste0(TAB,"DE_summary_shrunken.csv"), row.names=FALSE)

capture.output(sessionInfo(), file="results/diagnostics/R_sessionInfo_diagnostics.txt")
