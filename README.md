# scRNA-seq analysis — from quality control to a focused 933-cell subset

This repository is a visual, privacy-preserving record of a single-cell RNA-seq analysis performed in R with Seurat, without CellSweep. It presents the analytical path and its aggregate outputs, not a reusable data-release package.

> **Privacy notice.** No count matrix, sequencing file, serialized R object, sample identifier, condition label, gene name, reference object, or cell-level metadata is included. Every figure is an approved aggregate view with generic labels only.

## Result

The analysis converged on a focused subset of **933 cells** after quality review, reprocessing, and subset selection. The final representation retains the unsupervised structure required for downstream differential-expression screening, while withholding all study-specific identities.

| Focused 933-cell subset | Differential-expression screening |
|---|---|
| ![Focused 933-cell UMAP](assets/umap_final_subset_933.png) | ![Anonymized differential-expression volcano plots](assets/deg_volcano_plots.png) |

The UMAP groups transcriptionally similar cells in a two-dimensional representation. The volcano panels summarize the differential-expression screen by anonymized cluster. Neither plot names biological populations, samples, conditions, contrasts, or genes.

## Analysis sequence

### 1. Quality-control filtering

Technical QC metrics were reviewed across the pooled dataset, and cells failing the established quality criteria were removed before downstream reprocessing. The figure shows the aggregate change from **6,384** cells to **5,641** cells; no source-library or condition information is displayed.

![Quality-control filtering](assets/quality_control_filtering.png)

### 2. Marker dot plot

Cluster-associated features were reviewed alongside the QC results. The plot preserves the expression-pattern logic—colour for scaled expression and dot size for detection rate—while replacing every feature and cluster label with a generic identifier.

![Anonymized marker dot plot](assets/marker_dot_plot.png)

### 3. UMAP after the QC pass

The reprocessed, quality-filtered dataset was inspected in UMAP space to assess its unsupervised structure before narrowing to the focused population.

![UMAP after quality-control filtering](assets/umap_after_qc.png)

### 4. Focused subset of 933 cells

Cells selected for the focused analysis were reprocessed and embedded again. This is the final UMAP shown in the result section above; it is repeated here to keep the workflow chronological.

![Focused 933-cell UMAP](assets/umap_final_subset_933.png)

### 5. Differential-expression volcano plots

Differential-expression results were screened within the focused analysis. Each panel corresponds to an anonymized unsupervised cluster. Point identities and the compared conditions are intentionally withheld.

![Anonymized differential-expression volcano plots](assets/deg_volcano_plots.png)

## What is documented here

The complete annotated analysis narrative is available in [analysis/scRNAseq_workflow_sanitized.Rmd](analysis/scRNAseq_workflow_sanitized.Rmd). It records the quality assessment, dimensionality reduction, iterative filtering, marker review, subset analysis, and differential-expression stages without including restricted inputs or study-specific labels.

## Data protection and interpretation limits

- Restricted source data, per-cell tables, sample metadata, R objects, and arbitrary analysis outputs are excluded by [.gitignore](.gitignore).
- Labels in the figures are sanitized: clusters are generic and size-ranked, marker features are renumbered, and comparison labels are omitted.
- These aggregate plots document analytical capability and workflow progression. They are not a biological conclusion, a data release, or a clinical tool.
- Interpretation of the protected study remains confined to the authorized analysis environment.

## Main dependencies

R (>= 4.3), Seurat (v5), tidyverse, patchwork and, for nuclear-fraction calculation, DropletQC and Rsamtools.

## License

Add the license that applies before redistributing the code.
