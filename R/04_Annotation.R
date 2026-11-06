#!/usr/bin/env Rscript

# Remove excluded clusters and add the final cell-type annotations.

library(Seurat)

# 1. Input and output
combined_object <- readRDS(
  file.path("results", "rds", "marker_ready_object.rds")
)

annotation_plan <- read.csv(
  file.path("metadata", "final_cluster_annotation.csv")
)

output_dir <- file.path("results", "rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Keep clusters included in the core analysis
annotation_plan$cluster <- as.character(annotation_plan$cluster)
annotation_plan <- annotation_plan[
  annotation_plan$core_status == "INCLUDE",
]

cluster_column <- "cluster_harmony20_d10_r0_5"
cell_clusters <- as.character(combined_object[[cluster_column, drop = TRUE]])

cells_to_keep <- colnames(combined_object)[
  cell_clusters %in% annotation_plan$cluster
]

final_object <- subset(combined_object, cells = cells_to_keep)
final_clusters <- as.character(final_object[[cluster_column, drop = TRUE]])

# 3. Add major and subtype annotations
major_map <- setNames(
  annotation_plan$cell_type_major,
  annotation_plan$cluster
)

subtype_map <- setNames(
  annotation_plan$cell_type_subtype,
  annotation_plan$cluster
)

final_object$cell_type_major <- major_map[final_clusters]
final_object$cell_type_subtype <- subtype_map[final_clusters]

major_order <- c(
  "oRG", "gIPC", "nIPC", "Cycling gIPC",
  "OPC", "Cycling OPC", "NFOL", "EN", "IN",
  "Microglia", "VLMC", "Ependymal"
)

final_object$cell_type_major <- factor(
  final_object$cell_type_major,
  levels = major_order
)

Idents(final_object) <- final_object$cell_type_major

# 4. Save the final annotated object
saveRDS(
  final_object,
  file.path(output_dir, "final_annotated_object.rds")
)

