#!/usr/bin/env Rscript

# Create privacy-preserving figures from a local Seurat object.
# Usage: Rscript scripts/render_sanitized_figures.R /path/to/object.rds assets
# The input object and every intermediate output must remain outside version control.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: render_sanitized_figures.R <seurat_object.rds> <output_dir>")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

object_path <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(object_path)
required_columns <- c("nCount_RNA", "nFeature_RNA", "mitoRatio", "unintegrated_cluster_0.35")
if (!all(required_columns %in% colnames(seu[[]]))) {
  stop("The input object does not contain the expected QC and clustering metadata.")
}
if (!"umap.unintegrated" %in% names(seu@reductions)) {
  stop("The input object does not contain the expected UMAP reduction.")
}

# Generic labels prevent recovery of the study-specific cluster identity.
raw_cluster <- as.character(seu$unintegrated_cluster_0.35)
cluster_rank <- names(sort(table(raw_cluster), decreasing = TRUE))
generic_cluster <- factor(
  paste("Cluster", match(raw_cluster, cluster_rank)),
  levels = paste("Cluster", seq_along(cluster_rank))
)

emb <- as.data.frame(Embeddings(seu, reduction = "umap.unintegrated"))
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
emb$cluster <- generic_cluster

umap_plot <- ggplot(emb, aes(UMAP_1, UMAP_2, colour = cluster)) +
  geom_point(size = 0.12, alpha = 0.65) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(
    title = "Unsupervised cell clustering",
    subtitle = "Clusters are relabelled generically and ordered by size",
    colour = NULL,
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "sanitized_umap_clusters.png"), umap_plot,
       width = 8.5, height = 6.5, dpi = 240, bg = "white")

cluster_summary <- tibble(cluster = generic_cluster) |>
  count(cluster, name = "n_cells") |>
  mutate(cluster = forcats::fct_reorder(cluster, n_cells))

cluster_plot <- ggplot(cluster_summary, aes(cluster, n_cells)) +
  geom_col(fill = "#3264A8") +
  coord_flip() +
  labs(
    title = "Cells retained per unsupervised cluster",
    subtitle = "Generic cluster identifiers only",
    x = NULL,
    y = "Number of cells"
  ) +
  theme_classic(base_size = 11)

ggsave(file.path(output_dir, "sanitized_cluster_sizes.png"), cluster_plot,
       width = 7.5, height = 5.5, dpi = 240, bg = "white")

qc <- seu[[]] |>
  transmute(
    `UMI counts` = log10(nCount_RNA),
    `Detected features` = log10(nFeature_RNA),
    `Mitochondrial ratio` = mitoRatio
  ) |>
  tidyr::pivot_longer(everything(), names_to = "metric", values_to = "value")

qc_plot <- ggplot(qc, aes(metric, value, fill = metric)) +
  geom_violin(scale = "width", colour = NA, alpha = 0.9) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  labs(
    title = "Quality-control distributions after basic filtering",
    subtitle = "All cells pooled; no sample or condition identifiers are shown",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(output_dir, "sanitized_qc_distributions.png"), qc_plot,
       width = 11, height = 4.2, dpi = 240, bg = "white")
