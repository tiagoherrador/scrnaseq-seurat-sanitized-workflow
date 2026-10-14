#!/usr/bin/env Rscript

# Render privacy-reviewed, aggregate-only figures for the repository README.
# Inputs are local, access-controlled objects and must never be versioned.
# Usage:
# Rscript scripts/render_sanitized_figures.R \
#   <before_qc.rds> <after_qc.rds> <final_subset.rds> <degs.rds> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: render_sanitized_figures.R <before_qc.rds> <after_qc.rds> ",
    "<final_subset.rds> <degs.rds> <output_dir>"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

before_path <- args[[1]]
after_path <- args[[2]]
subset_path <- args[[3]]
degs_path <- args[[4]]
output_dir <- args[[5]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

before_qc <- readRDS(before_path)
after_qc <- readRDS(after_path)
final_subset <- readRDS(subset_path)
degs <- readRDS(degs_path)

required_qc <- c("nCount_RNA", "nFeature_RNA", "mitoRatio")
if (!all(required_qc %in% colnames(before_qc[[]])) ||
    !all(required_qc %in% colnames(after_qc[[]]))) {
  stop("Expected QC metadata are absent from one or more local objects.")
}

generic_clusters <- function(values) {
  ranked <- names(sort(table(as.character(values)), decreasing = TRUE))
  factor(
    paste("Cluster", match(as.character(values), ranked)),
    levels = paste("Cluster", seq_along(ranked))
  )
}

generic_umap <- function(seu, reduction, cluster_column, filename, title, subtitle) {
  if (!reduction %in% names(seu@reductions)) {
    stop("Expected dimensional reduction is absent from a local object.")
  }
  if (!cluster_column %in% colnames(seu[[]])) {
    stop("Expected clustering metadata are absent from a local object.")
  }

  coordinates <- as.data.frame(Embeddings(seu, reduction = reduction))
  colnames(coordinates)[1:2] <- c("UMAP_1", "UMAP_2")
  coordinates$cluster <- generic_clusters(seu[[cluster_column]][, 1])

  plot <- ggplot(coordinates, aes(UMAP_1, UMAP_2, colour = cluster)) +
    geom_point(size = 0.28, alpha = 0.72) +
    guides(colour = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
    labs(
      title = title,
      subtitle = subtitle,
      colour = NULL,
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, filename), plot,
         width = 9, height = 6.8, dpi = 240, bg = "white")
}

# 1. Quality-control filtering: pooled, aggregate views only.
qc_points <- bind_rows(
  before_qc[[]] |>
    transmute(
      stage = "Before QC filtering",
      library_size = log10(nCount_RNA),
      mitochondrial_ratio = mitoRatio
    ),
  after_qc[[]] |>
    transmute(
      stage = "After QC filtering",
      library_size = log10(nCount_RNA),
      mitochondrial_ratio = mitoRatio
    )
) |>
  group_by(stage) |>
  slice_sample(n = 4000L) |>
  ungroup() |>
  mutate(stage = factor(stage, levels = c("Before QC filtering", "After QC filtering")))

qc_scatter <- ggplot(qc_points, aes(library_size, mitochondrial_ratio)) +
  geom_point(colour = "#3264A8", alpha = 0.14, size = 0.42) +
  facet_wrap(~stage) +
  labs(
    title = "Quality-control filtering",
    subtitle = "Pooled technical metrics; no library, sample, or condition labels",
    x = "log10 total counts per cell",
    y = "Mitochondrial ratio"
  ) +
  theme_classic(base_size = 11)

qc_counts <- tibble(
  stage = factor(c("Before QC filtering", "After QC filtering"),
                 levels = c("Before QC filtering", "After QC filtering")),
  cells = c(ncol(before_qc), ncol(after_qc))
)

qc_funnel <- ggplot(qc_counts, aes(stage, cells, fill = stage)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(cells)), vjust = -0.35, size = 4) +
  labs(x = NULL, y = "Cells retained") +
  scale_fill_manual(values = c("#9CB9D9", "#3264A8")) +
  expand_limits(y = max(qc_counts$cells) * 1.12) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

ggsave(
  file.path(output_dir, "quality_control_filtering.png"),
  qc_scatter + qc_funnel + plot_layout(widths = c(2.1, 1)),
  width = 12.5, height = 5.2, dpi = 240, bg = "white"
)

# 2. Marker dot plot. Marker identities are selected locally then replaced by
# generic feature labels before plotting or writing any file.
marker_cluster <- "unintegrated_cluster_0.35"
after_qc <- JoinLayers(after_qc, assay = DefaultAssay(after_qc))
Idents(after_qc) <- after_qc[[marker_cluster]][, 1]
marker_table <- FindAllMarkers(
  after_qc,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  return.thresh = 0.05,
  verbose = FALSE
)

effect_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(marker_table))[[1]]
selected_markers <- marker_table |>
  group_by(cluster) |>
  arrange(desc(.data[[effect_column]]), .by_group = TRUE) |>
  slice_head(n = 2) |>
  ungroup() |>
  distinct(gene, .keep_all = TRUE) |>
  pull(gene)

if (length(selected_markers) < 2L) {
  stop("Insufficient local marker signal to create the aggregate dot plot.")
}

after_qc$privacy_cluster <- generic_clusters(after_qc[[marker_cluster]][, 1])
dot_plot <- DotPlot(
  after_qc,
  features = selected_markers,
  group.by = "privacy_cluster",
  assay = DefaultAssay(after_qc)
)
dot_plot$data$features.plot <- factor(
  paste("Feature", match(dot_plot$data$features.plot, selected_markers)),
  levels = paste("Feature", seq_along(selected_markers))
)
dot_plot$data$id <- factor(dot_plot$data$id, levels = levels(after_qc$privacy_cluster))

dot_plot <- dot_plot +
  scale_colour_gradient(low = "#E6EEF7", high = "#245A94") +
  scale_size(range = c(0.4, 6)) +
  labs(
    title = "Marker dot plot",
    subtitle = "Marker features and clusters are anonymized; dot size indicates detection rate",
    x = "Anonymized marker feature",
    y = "Unsupervised cluster",
    colour = "Scaled expression",
    size = "Detection rate"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "marker_dot_plot.png"), dot_plot,
       width = 11, height = 6.6, dpi = 240, bg = "white")

# 3. UMAP after the QC pass.
generic_umap(
  after_qc,
  reduction = "umap.unintegrated",
  cluster_column = "unintegrated_cluster_0.35",
  filename = "umap_after_qc.png",
  title = "UMAP after quality-control filtering",
  subtitle = "Unsupervised clusters are size-ranked and labelled generically"
)

# 4. Final focused subset (933 cells) after the downstream subset/reprocessing.
generic_umap(
  final_subset,
  reduction = "umap.harmony",
  cluster_column = "harmony_cluster_0.35",
  filename = "umap_final_subset_933.png",
  title = paste0("Focused subset after reprocessing (", ncol(final_subset), " cells)"),
  subtitle = "Final UMAP with anonymized, unsupervised cluster labels"
)

# 5. Differential-expression screening. Original feature names, contrasts and
# cluster codes are never exported; every visible label is generic.
if (!"0.35" %in% names(degs)) {
  stop("Expected DEG resolution is absent from the local DEG object.")
}
volcano_data <- degs[["0.35"]] |>
  filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Cluster)) |>
  mutate(
    cluster = generic_clusters(Cluster),
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    status = if_else(padj < 0.05 & abs(log2FoldChange) >= 0.25,
                     "Adjusted P < 0.05", "Not prioritized")
  )

volcano_plot <- ggplot(volcano_data, aes(log2FoldChange, neg_log10_padj, colour = status)) +
  geom_point(size = 0.34, alpha = 0.52) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  facet_wrap(~cluster, ncol = 3) +
  scale_colour_manual(values = c("Not prioritized" = "grey72", "Adjusted P < 0.05" = "#B2182B")) +
  labs(
    title = "Differential-expression screening",
    subtitle = "One panel per anonymized cluster; feature names and comparison labels are withheld",
    x = "Signed effect size",
    y = "−log10 adjusted P value",
    colour = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

volcano_height <- max(7, 2.45 * ceiling(nlevels(volcano_data$cluster) / 3) + 1.8)
ggsave(file.path(output_dir, "deg_volcano_plots.png"), volcano_plot,
       width = 12, height = volcano_height, dpi = 240, bg = "white")
