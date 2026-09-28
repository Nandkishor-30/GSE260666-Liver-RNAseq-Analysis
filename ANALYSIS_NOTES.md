# Detailed analysis notes

[Back to README](README.md)

These notes describe the supplied analysis and saved outputs. See [VALIDATION.md](VALIDATION.md) for verification limits.

## Quality Control

All 32 paired-end FASTQ files passed the main read-quality checks for per-base quality, per-sequence quality, adapter content, per-base N content, and sequence length. No substantial adapter contamination or poor-quality read tails were observed, so no trimming step was applied and STAR alignment was performed directly on the raw FASTQ files. (fastp appears in the initial software inventory but was not used in the final workflow; STAR soft-clips adapter-containing read ends during alignment.) High duplication levels were observed, which can occur in RNA-seq because highly expressed transcripts produce repeated reads. Reads were therefore not deduplicated.

## Reference and Alignment

The analysis used the GRCh38 primary assembly and GENCODE Release 48 annotation.

STAR genome indexing used:

`--sjdbOverhang 150` corresponding to 151 bp sequencing reads. Uniquely mapped read rates were approximately 91-95% across samples. Strandedness assessment indicated reverse-stranded libraries.

## Gene Quantification

featureCounts was run using paired-end fragment counting, reverse-stranded mode (`-s 2`), exon features, and `gene_id` grouping. Fragment assignment rates ranged from 69.7% to 81.4% per sample (mean 77.6%). Across all samples, 320.1M of 411.0M fragments were assigned to genes; the remainder were unassigned as multi-mapping (35.3M), ambiguous between overlapping features (19.5M), not overlapping any feature (23.7M), or unmapped (12.4M). Multi-mapping and ambiguous fragments were excluded under the counting settings used.

## Differential Expression

DESeq2 was used with:

`design = ~ condition` Low-expression genes were removed by requiring at least 10 counts in at least 4 samples.

- Genes before filtering: 78,894
- Genes after filtering: 17,679

Significance criteria:

- adjusted p-value < 0.05
- absolute log2 fold-change >= 1

| Comparison | DEGs | Up | Down |
|---|---:|---:|---:|
| NAFL vs Control | 17 | 10 | 7 |
| NASH vs Control | 474 | 184 | 290 |
| NASH vs NAFL | 28 | 10 | 18 |

Before applying the fold-change filter, the number of genes at adjusted p < 0.05 was 17 (NAFL vs Control), 825 (NASH vs Control) and 36 (NASH vs NAFL). The |log2FC| >= 1 requirement was applied in addition, since statistical significance alone can flag small effects of limited biological relevance.

## Exploratory Analysis

Variance-stabilizing transformation was used for visualization (`vst(blind = FALSE)`). The effect of using a blinded transform was not evaluated in this review.

PCA explained:

- PC1: 40.39%
- PC2: 10.56%
- PC1 + PC2: 50.95%

The first two components do not cleanly separate the three conditions: Control, NAFL and NASH samples are interspersed along PC1, and the sample-to-sample distance heatmap likewise does not cluster strictly by condition. This is expected for human liver tissue, where inter-individual variation is large relative to the disease effect, and it is consistent with the modest DEG counts for NAFL vs Control. NASH samples nonetheless occupy a partially distinct region, and the supervised DESeq2 analysis, which models condition explicitly rather than relying on unsupervised structure, identifies condition-associated differences that require independent validation.

## Functional Enrichment

Enrichment uses genes mapped to Entrez identifiers through `org.Hs.eg.db`. Unmapped genes are excluded and duplicate Entrez identifiers are collapsed for GSEA. Mapping counts are printed during execution; these counts can differ from the number of input Ensembl identifiers because mappings are not necessarily one-to-one. GSEA uses a fixed seed (`set.seed(123)`, `seed = TRUE`). Exact reproducibility across software versions or parallel backends has not been established. Over-representation analysis was run separately for upregulated and downregulated genes. For NASH vs Control, 413 of 474 DEGs mapped to Entrez identifiers (170 up, 243 down). Upregulated genes yielded 17 significantly enriched GO BP terms.

Downregulated genes, and the combined DEG list, yielded none. No significant terms were detected for the downregulated or combined lists under these settings. This can reflect limited power, gene-set coverage, heterogeneity, or other factors; it does not establish a single explanation. GO Biological Process enrichment was performed using clusterProfiler and `org.Hs.eg.db`.

For NASH vs Control, upregulated genes showed enrichment for processes including:

- ribonucleoprotein complex biogenesis
- ribosome biogenesis
- rRNA processing
- ribosome assembly
- leukocyte aggregation
- antimicrobial humoral responses
- regulation of apoptotic signaling

GSEA identified:

| Comparison | Significant GO-BP GSEA terms |
|---|---:|
| NAFL vs Control | 307 |
| NASH vs Control | 660 |
| NASH vs NAFL | 647 |

GO terms are hierarchical and redundant, so these were interpreted as broader biological themes rather than independent pathways.

Major recurring themes included:

- ribosome biogenesis and translation
- RNA and rRNA processing
- amino-acid and organic-acid metabolism
- cellular and ion homeostasis
- immune and antimicrobial responses
- apoptotic and stress signaling
- extracellular matrix and tissue organization

## Biological Interpretation

NAFL showed relatively modest differential expression at the individual-gene level (17 genes past the combined thresholds), but GSEA revealed a strong and internally consistent signal: the most significantly depleted gene sets in NAFL relative to Control were all amino-acid catabolic processes (alpha-amino acid catabolism NES -2.58, amino acid catabolism NES -2.47, L-amino acid catabolism NES -2.39). This contrast between a near-null per-gene result and a coherent pathway-level result illustrates why GSEA complements, rather than duplicates, individual-gene testing.

NASH showed substantially greater transcriptional alteration, with strong enrichment of ribosome biogenesis, RNA processing, translation-related programs, and immune/stress-associated processes. Metabolic processes including amino-acid and organic-acid catabolism showed reduced enrichment relative to healthy Control samples. The NASH vs NAFL comparison additionally highlighted changes in RNA-processing, translational, extracellular-matrix, and tissue-organization programs.

## Pipeline Validation

Three exploratory diagnostics were examined. These checks do not rule out confounding or technical artefacts.

### Library composition

Script `11` estimates biotype fractions using raw counts for the 17,679 genes retained after expression filtering. These fractions describe counts within that filtered gene set, not all sequenced reads or all assigned fragments. The reported protein-coding fraction was approximately 96.96 per cent (mean across samples); retained rRNA-derived counts were predominantly mitochondrial. This analysis alone cannot establish the efficiency of rRNA depletion or exclude technical explanations for enrichment results.

### Model diagnostics

Dispersion estimates shrank toward the fitted trend as expected. P-value distributions were approximately uniform with an excess near zero, but this visual check alone does not establish test calibration. Log2 fold-change shrinkage was applied using apeglm for the two contrasts that correspond to model coefficients (NAFL vs Control, NASH vs Control) and ashr for NASH vs NAFL, which is not a model coefficient.

Shrinkage materially changed the number of genes passing the combined significance and fold-change criteria:

| Comparison | Genes passing, unshrunken | Genes passing, shrunken |
|---|---:|---:|
| NAFL vs Control | 17 | 15 |
| NASH vs Control | 474 | 337 |
| NASH vs NAFL | 28 | 1 |

Shrinkage reduces the number of genes exceeding the additional fold-change threshold, especially for NASH vs NAFL. This highlights uncertainty in effect-size estimates in this small cohort. The affected genes were not individually characterized here. The diagnostic script retains the original DESeq2 p-values and adjusted p-values and reports shrunken effect sizes separately. The Wald statistic uses the estimated coefficient divided by its standard error. Shrinkage does not replace the original hypothesis tests in these tables. Both estimates are provided in `results/diagnostics/`.

### Marker genes

Expression of established NAFLD/NASH markers was examined directly. Mean expression increased from Control through NAFL to NASH for AKR1B10 (5.30, 6.00, 8.07 VST), SPP1 (9.34, 9.49, 10.21), CCL2 (7.59, 8.09, 9.04) while KRT19 was non-monotonic (6.26, 5.84, 6.92). Fibrosis and extracellular-matrix markers did not follow this pattern. COL1A1, COL1A2, COL3A1 and TIMP1 peaked in NAFL rather than NASH. One possible explanation is that NASH is diagnosed on steatosis, inflammation and hepatocyte ballooning, whereas fibrosis is staged on a separate axis; a NASH cohort can therefore carry little fibrosis.

With four NASH samples and no fibrosis staging in the available metadata, this cannot be resolved here, and the discordance is reported rather than smoothed over. This is not an isolated marker-level observation. GSEA of the same NASH versus NAFL expression data shows coordinated suppression of collagen fibril organization, collagen metabolic process and extracellular matrix organization, and the same extracellular matrix terms appear among the suppressed sets in NASH versus Control. Marker-gene and pathway-level evidence therefore agree, but neither establishes fibrosis stage nor rules out technical artefacts or confounding.

### GC-content review — 27 September 2026

The saved MultiQC report records 9 passes, 19 warnings and 4 failures for per-sequence GC content across 32 FASTQ files. Failures occurred in SRR28195703_2 and SRR28195715_2 (Control), and SRR28195707_2 and SRR28195708_2 (NAFL). Inspection of the embedded distributions showed shoulders/multiple peaks around 45–56% GC. The flagged read 2 distributions broadly resembled their read 1 partners and the cohort distributions. This observation does not establish the cause or exclude contamination or library-preparation bias. No samples were excluded and no additional trimming was performed based solely on these flags. Their cause remains unresolved.

