# Setup and execution

Run all commands from the main project folder containing `README.md`, `environment.yml`, and `scripts/`. This folder is called the repository root.

For a shorter rerun, start with the included count matrix using “Reanalyse the supplied counts.” To repeat downloading, quality control, alignment, and counting, follow “Rebuild from raw reads” after installing the dependencies.

## Upload to GitHub

1. Extract the ZIP on your computer.
2. Create an empty GitHub repository.
3. Choose **Add file → Upload files** and upload the extracted files and folders.
4. Ensure README.md, SETUP.md, environment.yml, scripts/, counts/, metadata/, results/, figures/, and qc/ are at the repository root. Include `.gitignore` and `.gitattributes` (hidden files may need to be shown in your file manager).
5. Commit the files. Do not upload only the ZIP: GitHub does not unpack it.

No Git LFS is needed for the included files. GitHub displays PNG figures; open the downloaded MultiQC HTML locally to use its interactive report.

## Dependencies

Use Linux or a Linux HPC environment for the raw-read pipeline. A human STAR index requires substantial memory; allocate appropriate compute and disk space before downloading all runs. Thread count is configurable through `THREADS`.

With Conda installed, from the repository root:

```bash
conda env create -f environment.yml
conda activate gse260666
Rscript scripts/install_R_packages.R
```

The installation recipe targets R 4.4 and Bioconductor 3.20. The counts-based validation used R 4.4.3 and ggplot2 3.5.2 after resolving dependency issues. The revised installation recipe has not yet been independently tested from scratch and is not an exact environment lockfile.

Original software versions are recorded in `metadata/R_sessionInfo_final.txt` and `metadata/software_versions_initial.txt`. Records from the validation environment are available in `validation/`. Both `apeglm` and `ashr` are required for the shrinkage diagnostics. Installation may require compilers and system libraries.

On HPC, existing modules can instead provide the tools. Load them before running scripts. Scripts use executables on PATH; `STAR_BIN` and `FEATURECOUNTS_BIN` can override the STAR and featureCounts executables.

## Reanalyse the supplied counts

```bash
conda activate gse260666
bash scripts/run_from_counts.sh
```

This runs scripts 06–14, including generation of the two RDS objects used by later steps. It needs no FASTQ, BAM, or reference download. All R scripts use paths relative to the repository root. New results overwrite saved outputs; work in a separate checkout when comparing results.

## Rebuild from raw reads

```bash
export THREADS=8
bash scripts/00_prepare_reference.sh
bash scripts/01_download_fastq.sh
bash scripts/02_fastqc_multiqc.sh
```

Inspect `qc/multiqc_report.html` before alignment. The original analysis did not trim reads; reassess this decision if your input data differ.

```bash
bash scripts/03_reference_STAR_index.sh
bash scripts/04_STAR_alignment.sh
python3 scripts/summarize_STAR_QC.py
```

Inspect `results/STAR_alignment_QC.tsv` and `results/STAR_strandedness_check.tsv`. Script 05 uses reverse-stranded counting (`-s 2`) as in the supplied analysis. Confirm this remains appropriate before continuing.

```bash
bash scripts/05_featureCounts.sh
bash scripts/run_from_counts.sh
```

The reference preparation script downloads the release-48 primary-assembly FASTA and comprehensive GTF from GENCODE. The annotation helper extracts gene ID, symbol and biotype from GTF gene records, preserving versioned IDs. If the supplied annotation disagrees, it stops for investigation rather than replacing it silently. The STAR helper regenerates the two summary tables from STAR logs and ReadsPerGene files; those raw logs are not included in this compact package. These helpers were added during packaging, not recovered from the original historical run.

## Restart behavior

The downloader validates gzip integrity for existing pairs and uses a temporary staging directory for new conversions. Failed staging directories may remain under `tmp/` for investigation. STAR alignment is skipped only when a nonempty BAM and completion marker are both present; a partial BAM alone is not accepted. No script intentionally deletes raw reads or aligned BAMs.

## September 2026 review updates

The DESeq2 rerun matched all six original result tables within numerical tolerance of 1e-8 after sorting by gene ID. PCA labels and condition colors were improved.

The installer restores ggplot2 3.5.2 after dependency installation. It passed in the repaired environment; a fresh installation of the revised installer remains untested.

Optional sensitivity analysis: run `Rscript scripts/GSEA_outlier_sensitivity.R` from the repository root after the counts-based workflow.

Excluding outlier-filtered genes from GSEA retained 91.2–96.8% of the original significant terms, with no direction reversals among shared terms. Leading themes persisted; individual term lists changed. Software and sampling differences may also contribute.

See [VALIDATION.md](VALIDATION.md) and [the comparison table](validation/GSEA_sensitivity_summary.csv). GC-content flags remain unresolved and are documented in [ANALYSIS_NOTES.md](ANALYSIS_NOTES.md).
