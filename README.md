# Reproducible scRNA-seq workflow (data excluded)

This repository contains an anonymized version of an original single-cell RNA-seq analysis conducted with R and Seurat, without CellSweep. It provides the analysis notebook and an example configuration; it contains no data, sample identifiers, local paths, results, or serialized objects.

## Scope

The workflow covers 10x matrix import, quality control, optional nuclear-fraction estimation, normalization, PCA, UMAP, clustering, marker detection, iterative filtering, reference-based annotation, and subsetting. Sample, group, and reference names have been replaced with generic placeholders.

## Usage

1. Create a local data/ directory, which is not versioned, and place matrices and any required BAM files there.
2. Copy config/example_samples.csv to config/samples.csv and fill in approved local paths and identifiers for your environment.
3. Review the QC parameters and optional sections of the notebook.
4. Run or render analysis/scRNAseq_workflow_sanitized.Rmd locally.

The notebook uses eval = FALSE to prevent accidental execution on private data. Deliberately enable only the chunks needed in your local copy.

## Data protection

- Never add matrices, BAM/FASTQ files, metadata, serialized R objects, per-cell results, or figures derived from restricted data.
- Keep local paths and access credentials out of version control.
- Review git status before every git add.

## Main dependencies

R (>= 4.3), Seurat (v5), tidyverse, patchwork and, for nuclear-fraction calculation, DropletQC and Rsamtools.

## License

Add the license that applies before redistributing the code.

