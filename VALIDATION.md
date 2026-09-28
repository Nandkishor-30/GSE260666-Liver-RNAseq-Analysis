# Validation and reproducibility

## Counts-based rerun — 27 September 2026

### featureCounts input checks — 27 September 2026

Script 05 now selects BAMs from metadata/SRR_accessions.txt rather than a directory wildcard. It rejects invalid or duplicate SRR accessions, an empty accession list, missing or empty BAMs, and missing alignment completion markers. Shell syntax validation passed. These checks have not been tested with real BAM inputs; markers do not prove BAM integrity.

The counts-based workflow (scripts 06–13) was rerun on the university Linux server in a newly created Conda environment, after resolving installation dependencies and using ggplot2 3.5.2 with R 4.4.3 and Bioconductor 3.20.

### Confirmed results

- All eight required R packages loaded successfully.
- Sample matching confirmed 6 Control, 6 NAFL and 4 NASH samples.
- Filtering retained 17,679 of 78,894 genes.
- All six DESeq2 tables (all-gene and significant-DEG tables for three comparisons) matched the original tables using R all.equal with tolerance = 1e-8, after sorting by gene ID.
- All six enrichment CSV tables matched using the same tolerance, after sorting by term ID where available.
- DEG totals were 17, 474 and 28.
- GSEA term totals were 307, 660 and 647.
- Reported shrunken-threshold totals were 15, 337 and 1.
- The run log reached the final marker-panel script without a visible error. The shell exit status was not captured.

### Scope and remaining limitations

This validates the 12 compared result tables from the included count matrix. Raw-read downloading, QC, alignment and counting were not rerun. Other output files were not comprehensively compared.

Warnings included omitted plotting rows, unmapped gene identifiers and tied GSEA ranking statistics. Their implications still require review.

The setup recipe now incorporates the dependency fixes used in this environment. The revised recipe has not itself been independently tested from scratch, and is not an exact lockfile.

### Included validation records

The validation folder contains:

-  counts_validation_run.log: counts-based analysis run log.
-  DEG_comparison.log: comparison of the regenerated DESeq2 tables with the original tables.
-  enrichment_comparison.log: comparison of the regenerated enrichment tables with the original tables.
-  conda_explicit_20260927.txt: recorded Conda package specifications.
-  R_packages_20260927.csv: installed R package inventory.
-  R_sessionInfo_20260927.txt: recorded R session information.
-  GSEA_outlier_sensitivity.log, GSEA_sensitivity_summary.csv, and GSEA_sensitivity_sessionInfo.txt: sensitivity-analysis records.

The original-output backup mentioned in the local validation history is not included in this project package.

## Packaging review — 26 September 2026

### Verified against supplied files

- ZIP integrity and count-matrix gzip integrity passed.
- All 16 count samples match supplied sample metadata and SRA run accessions.
- All SRR-to-GSM mappings agree with the supplied SRA table.
- Counts are nonnegative integers: 78,894 genes by 16 samples.
- The expression filter retains 17,679 genes.
- Count-column totals equal featureCounts assigned-fragment totals.
- Assignment percentages range from 69.7466% to 81.4252% (mean 77.6170%).
- DEG totals, directions, and significant-gene identifiers agree with all-gene tables:

| Comparison | DEGs | Up | Down | padj < 0.05 | Shrunken threshold | GSEA terms |
|---|---:|---:|---:|---:|---:|---:|
| NAFL vs Control | 17 | 10 | 7 | 17 | 15 | 307 |
| NASH vs Control | 474 | 184 | 290 | 825 | 337 | 660 |
| NASH vs NAFL | 28 | 10 | 18 | 36 | 1 | 647 |

- Gene annotation identifiers are unique and cover all count-matrix genes.
- All distributed shell scripts pass `bash -n`; added Python helpers compile.
- No common GitHub token, AWS access-key, or private-key patterns were found in the inspected text files. This is a limited scan, not a security certification.
- No included file exceeds GitHub's browser-upload size limit.

### Changes in this package

Portable project paths and executable lookup; explicit dependencies; reference and annotation preparation helpers; STAR QC extraction helper; counts-only runner; safer interrupted-download and alignment handling; explicit shrinkage failures; stricter sample checks; fixed plotting seed for future marker plots; corrected README explanations and bounded biological interpretations; setup instructions; ignore rules for environments and intermediate R objects.

Original saved results, figures, count data and metadata were preserved byte for byte. New helper scripts do not establish the provenance of historical outputs. SHA256SUMS covers the distributed files other than itself.

### Not verified during the 26 September review

R and the bioinformatics executables were unavailable in the review environment. No dependency environment was solved, and no raw-data pipeline, DESeq2, shrinkage, or enrichment analysis was rerun. R changes received source review, not runtime validation. GENCODE download URLs were checked against its release page; the large references were not downloaded. Live GEO was blocked, so condition labels were not independently confirmed against the original records. The original MultiQC claims and every biological interpretation were not independently established. Neither numerical consistency nor the presence of diagnostics proves absence of confounding or technical artefacts.

### Documentation sources

- https://www.gencodegenes.org/human/release_48.html
- https://bioconductor.org/news/bioc_3_20_release/
- https://bioconductor.org/packages/DESeq2/
- https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github

## GSEA outlier sensitivity — 27 September 2026

Diagnostic comparison confirmed that 24 genes per contrast had p-values withheld by DESeq2 Cook's-distance filtering. The sensitivity run excluded genes with missing raw p-values from GSEA rankings, retaining genes excluded only by independent filtering. Ranked Entrez identifiers decreased from 15,583 to 15,564.

Relative to the saved original results, 91.2%, 94.1% and 96.8% of significant terms were retained for NAFL vs Control, NASH vs Control and NASH vs NAFL, respectively. No shared terms reversed enrichment direction. Leading metabolic and ribosome/RNA-processing themes persisted.

This supports stability of the leading themes, not every individual term. Differences may also reflect software versions and GSEA sampling variability because the comparison used the saved original run rather than a contemporaneous baseline.

Summary: validation/GSEA_sensitivity_summary.csv Results: results/enrichment_outlier_sensitivity/ Run log: validation/GSEA_outlier_sensitivity.log
