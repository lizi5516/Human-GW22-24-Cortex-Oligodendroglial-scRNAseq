#!/usr/bin/env Rscript

# Identify positive marker genes for every cluster.

library(Seurat)

set.seed(123)

# 1. Input and output
combined_object <- readRDS(
  file.path("results", "rds", "integrated_object.rds")
)

rds_dir <- file.path("results", "rds")
table_dir <- file.path("results", "tables")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Log-normalize the RNA assay
DefaultAssay(combined_object) <- "RNA"
combined_object <- JoinLayers(combined_object, assay = "RNA")
combined_object <- NormalizeData(combined_object, scale.factor = 10000)

Idents(combined_object) <- combined_object$cluster_harmony20_d10_r0_5

# 3. Find positive markers
cluster_markers <- FindAllMarkers(
  combined_object,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  max.cells.per.ident = 5000,
  random.seed = 123
)

# 4. Save the object and marker table
saveRDS(
  combined_object,
  file.path(rds_dir, "marker_ready_object.rds")
)

write.csv(
  cluster_markers,
  file.path(table_dir, "cluster_markers.csv"),
  row.names = FALSE
)

