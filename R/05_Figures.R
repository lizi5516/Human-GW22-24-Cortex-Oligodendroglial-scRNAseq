#!/usr/bin/env Rscript

# Draw the final annotation UMAP and marker heatmap.

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# 1. Input and output
final_object <- readRDS(
  file.path("results", "rds", "final_annotated_object.rds")
)

marker_table <- read.csv(
  file.path("metadata", "marker_heatmap_genes.csv")
)

output_dir <- file.path("results", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Cell-type order and colors
major_order <- c(
  "oRG", "gIPC", "nIPC", "Cycling gIPC",
  "OPC", "Cycling OPC", "NFOL", "EN", "IN",
  "Microglia", "VLMC", "Ependymal"
)

major_colors <- c(
  "oRG" = "#7C4DFF",
  "gIPC" = "#AB47BC",
  "nIPC" = "#EC407A",
  "Cycling gIPC" = "#F06292",
  "OPC" = "#1976D2",
  "Cycling OPC" = "#42A5F5",
  "NFOL" = "#26C6DA",
  "EN" = "#EF5350",
  "IN" = "#66BB6A",
  "Microglia" = "#8D6E63",
  "VLMC" = "#FFB300",
  "Ependymal" = "#78909C"
)

final_object$cell_type_major <- factor(
  final_object$cell_type_major,
  levels = major_order
)

# 3. Final annotation UMAP
umap_plot <- DimPlot(
  final_object,
  reduction = "umap_harmony20_d10",
  group.by = "cell_type_major",
  cols = major_colors,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15
) +
  labs(title = "Human OPC single-cell atlas") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )

ggsave(
  file.path(output_dir, "final_annotation_umap.pdf"),
  umap_plot,
  width = 10,
  height = 8
)

ggsave(
  file.path(output_dir, "final_annotation_umap.png"),
  umap_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# 4. Average marker expression by major cell type
marker_genes <- unique(marker_table$gene)
marker_genes <- marker_genes[marker_genes %in% rownames(final_object)]

DefaultAssay(final_object) <- "RNA"
Idents(final_object) <- final_object$cell_type_major

average_expression <- AverageExpression(
  final_object,
  features = marker_genes,
  group.by = "cell_type_major",
  layer = "data"
)$RNA

average_expression <- average_expression[
  marker_genes,
  major_order,
  drop = FALSE
]

# Convert each marker gene to a row Z-score.
z_score_matrix <- t(scale(t(average_expression)))
z_score_matrix[is.na(z_score_matrix)] <- 0
z_score_matrix <- pmax(pmin(z_score_matrix, 2), -2)

heatmap_data <- as.data.frame(z_score_matrix) |>
  rownames_to_column("gene") |>
  pivot_longer(
    cols = -gene,
    names_to = "cell_type_major",
    values_to = "z_score"
  ) |>
  mutate(
    gene = factor(gene, levels = rev(marker_genes)),
    cell_type_major = factor(cell_type_major, levels = major_order)
  )

# 5. Marker heatmap
heatmap_plot <- ggplot(
  heatmap_data,
  aes(x = cell_type_major, y = gene, fill = z_score)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    limits = c(-2, 2),
    name = "Row Z-score"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic")
  )

ggsave(
  file.path(output_dir, "marker_heatmap.pdf"),
  heatmap_plot,
  width = 10,
  height = 12
)

ggsave(
  file.path(output_dir, "marker_heatmap.png"),
  heatmap_plot,
  width = 10,
  height = 12,
  dpi = 300
)

