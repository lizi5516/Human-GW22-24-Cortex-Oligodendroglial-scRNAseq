#!/usr/bin/env Rscript

# Read Cell Ranger matrices and perform basic QC for all samples.

library(Seurat)

set.seed(123)

# 1. Input and output
cellranger_dir <- Sys.getenv(
  "SC_CELLRANGER_ROOT",
  unset = file.path("data", "cellranger")
)

output_dir <- file.path("results", "rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_metadata <- read.csv(
  file.path("metadata", "sample_metadata.csv")
)

# 2. QC thresholds
min_features <- 700
min_counts <- 1000
max_percent_mt <- 15

# 3. Read and filter the 16 samples
sample_list <- lapply(seq_len(nrow(sample_metadata)), function(i) {
  sample_id <- sample_metadata$sample_id[i]
  matrix_dir <- file.path(
    cellranger_dir,
    sample_metadata$matrix_relative_path[i]
  )

  counts <- Read10X(matrix_dir)

  # Keep only the gene-expression matrix if Read10X returns multiple assays.
  if (is.list(counts)) {
    counts <- counts[["Gene Expression"]]
  }

  sample_object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )

  sample_object <- RenameCells(
    sample_object,
    add.cell.id = sample_id
  )

  sample_object$sample_id <- sample_id
  sample_object$donor_id <- sample_metadata$donor_id[i]
  sample_object$sort_marker <- sample_metadata$sort_marker[i]
  sample_object$gestational_week <- sample_metadata$gestational_week[i]
  sample_object$batch <- sample_metadata$batch[i]

  sample_object[["percent.mt"]] <- PercentageFeatureSet(
    sample_object,
    pattern = "^MT-"
  )

  subset(
    sample_object,
    subset = nFeature_RNA >= min_features &
      nCount_RNA >= min_counts &
      percent.mt <= max_percent_mt
  )
})

names(sample_list) <- sample_metadata$sample_id

# 4. Save the QC-filtered samples
saveRDS(
  sample_list,
  file.path(output_dir, "qc_filtered_samples.rds")
)

