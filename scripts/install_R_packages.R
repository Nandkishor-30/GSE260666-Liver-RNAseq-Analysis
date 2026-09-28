#!/usr/bin/env Rscript

if (getRversion() < "4.4" || getRversion() >= "4.5")
    stop("Use R 4.4.x with Bioconductor 3.20")

# Each phase runs in a fresh R process to avoid stale loaded namespaces.
run_phase <- function(code) {
    script <- tempfile(fileext = ".R")
    on.exit(unlink(script))
    writeLines(code, script)
    status <- system2(
        file.path(R.home("bin"), "Rscript"),
        shQuote(script)
    )
    if (status != 0L)
        stop("Installation phase failed. Review the output above.")
}

run_phase('
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("mgcv", quietly = TRUE))
    stop("Install r-mgcv through the Conda environment first")

BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)

# Bootstrap dependencies. ggtree may need the compatibility repair below.
BiocManager::install(
    c("DESeq2", "clusterProfiler", "org.Hs.eg.db",
      "enrichplot", "apeglm"),
    ask = FALSE, update = FALSE
)

for (p in c("pheatmap", "ashr", "ggrepel")) {
    if (!p %in% rownames(installed.packages()))
        install.packages(p)
}
')

run_phase('
target <- "3.5.2"
installed <- installed.packages()
if (!"ggplot2" %in% rownames(installed) ||
    installed["ggplot2", "Version"] != target) {
    install.packages(
        paste0(
            "https://cran.r-project.org/src/contrib/Archive/ggplot2/",
            "ggplot2_", target, ".tar.gz"
        ),
        repos = NULL, type = "source"
    )
}
stopifnot(as.character(packageVersion("ggplot2")) == target)
')

run_phase('
BiocManager::install(
    c("ggtree", "enrichplot", "clusterProfiler"),
    ask = FALSE, update = FALSE, dependencies = FALSE
)
')

run_phase('
required <- c(
    "DESeq2", "ggplot2", "pheatmap", "clusterProfiler",
    "org.Hs.eg.db", "enrichplot", "apeglm", "ashr",
    "ggtree", "ggrepel"
)
missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing))
    stop("Missing or unloadable packages: ", paste(missing, collapse = ", "))

stopifnot(as.character(packageVersion("ggplot2")) == "3.5.2")
cat("All required packages load successfully.\\n")
')
