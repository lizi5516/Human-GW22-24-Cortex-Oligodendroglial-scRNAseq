#!/usr/bin/env Rscript

# Merge samples and perform SCTransform, PCA, Harmony, clustering and UMAP.

library(Seurat)
library(harmony)

set.seed(123)

# 1. Input and output
sample_list <- readRDS(
  file.path("results", "rds", "qc_filtered_samples.rds")
)

output_dir <- file.path("results", "rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Merge all samples
combined_object <- merge(
  x = sample_list[[1]],
  y = sample_list[-1],
  project = "Human_OPC"
)

# 3. Select 3,000 variable genes
DefaultAssay(combined_object) <- "RNA"

combined_object <- FindVariableFeatures(
  combined_object,
  selection.method = "vst",
  nfeatures = 3000
)

variable_genes <- VariableFeatures(combined_object)

# 4. Run SCTransform v2
combined_object <- SCTransform(
  combined_object,
  assay = "RNA",
  new.assay.name = "SCT",
  vst.flavor = "v2",
  residual.features = variable_genes,
  return.only.var.genes = TRUE,
  conserve.memory = TRUE,
  do.correct.umi = FALSE,
  seed.use = 123
)

# 5. Run PCA
combined_object <- RunPCA(
  combined_object,
  assay = "SCT",
  npcs = 50,
  seed.use = 123
)

# 6. Correct donor effects with Harmony
combined_object <- RunHarmony(
  combined_object,
  group.by.vars = "donor_id",
  reduction.use = "pca",
  dims.use = 1:20,
  reduction.save = "harmony_pca20"
)

# 7. Build the graph and identify clusters
combined_object <- FindNeighbors(
  combined_object,
  reduction = "harmony_pca20",
  dims = 1:10,
  k.param = 30,
  prune.SNN = 1 / 15,
  nn.method = "annoy",
  graph.name = c("harmony20_d10_nn", "harmony20_d10_snn")
)

combined_object <- FindClusters(
  combined_object,
  graph.name = "harmony20_d10_snn",
  resolution = 0.5,
  algorithm = 4,
  random.seed = 123,
  cluster.name = "cluster_harmony20_d10_r0_5"
)

# 8. Run UMAP
combined_object <- RunUMAP(
  combined_object,
  reduction = "harmony_pca20",
  dims = 1:10,
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  metric = "cosine",
  seed.use = 123,
  reduction.name = "umap_harmony20_d10",
  reduction.key = "UMAP_"
)

Idents(combined_object) <- combined_object$cluster_harmony20_d10_r0_5

# 9. Save the integrated object
saveRDS(
  combined_object,
  file.path(output_dir, "integrated_object.rds")
)

