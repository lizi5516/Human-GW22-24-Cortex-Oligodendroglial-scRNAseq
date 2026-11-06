#!/usr/bin/env Rscript

# Calculate and plot cell-type composition for the enrichment samples.

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggrastr)

# 1. Input and output
final_object <- readRDS(
  file.path("results", "rds", "final_annotated_object.rds")
)

sample_mapping <- read.csv(
  file.path("metadata", "enrichment_sample_mapping.csv")
)

figure_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Extract metadata and UMAP coordinates
cell_metadata <- final_object[[]] |>
  rownames_to_column("cell_id") |>
  left_join(sample_mapping, by = "sample_id")

umap_table <- as.data.frame(
  Embeddings(final_object, reduction = "umap_harmony20_d10")
)
colnames(umap_table)[1:2] <- c("UMAP_1", "UMAP_2")
umap_table <- rownames_to_column(umap_table, "cell_id")

cell_metadata <- left_join(
  cell_metadata,
  umap_table,
  by = "cell_id"
)

treatment_order <- c("Unenriched", "PDGFRA", "EGFR", "LIFR")

major_order <- c(
  "oRG", "gIPC", "nIPC", "Cycling gIPC",
  "OPC", "Cycling OPC", "NFOL", "EN", "IN",
  "Microglia", "VLMC", "Ependymal"
)

cell_metadata <- cell_metadata |>
  mutate(
    treatment = factor(treatment, levels = treatment_order),
    cell_type_major = factor(cell_type_major, levels = major_order)
  )

# 3. Calculate cell-type composition
composition_table <- cell_metadata |>
  count(
    donor_label,
    treatment,
    cell_type_major,
    name = "n_cells"
  ) |>
  group_by(donor_label, treatment) |>
  complete(
    cell_type_major = major_order,
    fill = list(n_cells = 0)
  ) |>
  mutate(
    percentage = 100 * n_cells / sum(n_cells)
  ) |>
  ungroup()

write.csv(
  composition_table,
  file.path(table_dir, "enrichment_cell_type_composition.csv"),
  row.names = FALSE
)

# 4. Cell-type colors
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

# 5. Donor cell-composition plot
composition_plot <- ggplot(
  composition_table,
  aes(x = treatment, y = percentage, fill = cell_type_major)
) +
  geom_col() +
  facet_wrap(~ donor_label, scales = "free_x") +
  scale_fill_manual(values = major_colors) +
  labs(x = NULL, y = "Cell composition (%)", fill = "Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(figure_dir, "donor_cell_composition.pdf"),
  composition_plot,
  width = 12,
  height = 5
)

ggsave(
  file.path(figure_dir, "donor_cell_composition.png"),
  composition_plot,
  width = 12,
  height = 5,
  dpi = 300
)

# 6. Paired comparison for Donor 1 and Donor 2
paired_table <- composition_table |>
  filter(donor_label %in% c("Donor 1", "Donor 2"))

paired_plot <- ggplot(
  paired_table,
  aes(
    x = treatment,
    y = percentage,
    group = donor_label,
    color = donor_label
  )
) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ cell_type_major, scales = "free_y", ncol = 4) +
  scale_color_manual(
    values = c("Donor 1" = "#1B9E77", "Donor 2" = "#D95F02")
  ) +
  labs(x = NULL, y = "Cell composition (%)", color = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(figure_dir, "paired_donor_comparison.pdf"),
  paired_plot,
  width = 12,
  height = 9
)

ggsave(
  file.path(figure_dir, "paired_donor_comparison.png"),
  paired_plot,
  width = 12,
  height = 9,
  dpi = 300
)

# 7. UMAP split by enrichment treatment
umap_plot <- ggplot(
  cell_metadata,
  aes(x = UMAP_1, y = UMAP_2, color = treatment)
) +
  geom_point_rast(size = 0.15) +
  facet_wrap(~ treatment, nrow = 1) +
  scale_color_manual(
    values = c(
      "Unenriched" = "#666666",
      "PDGFRA" = "#377EB8",
      "EGFR" = "#E41A1C",
      "LIFR" = "#984EA3"
    )
  ) +
  coord_equal() +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  file.path(figure_dir, "umap_by_treatment.pdf"),
  umap_plot,
  width = 14,
  height = 4
)

ggsave(
  file.path(figure_dir, "umap_by_treatment.png"),
  umap_plot,
  width = 14,
  height = 4,
  dpi = 300
)

