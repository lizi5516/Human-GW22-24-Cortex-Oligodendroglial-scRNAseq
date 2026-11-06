# Human-GW22-24-Cortex-Oligodendroglial-scRNAseq

Code and analysis workflows for a gliogenic-lineage enriched single-cell transcriptomic dataset of the developing human cortex.



## Overview

This repository contains the analysis code accompanying the dataset:

**A Gliogenic-lineage Enriched Single-cell Transcriptomic Dataset of Developing Human Cortex**

The dataset comprises single-cell RNA-sequencing data from the human gestational week 22–24 cortex. The analysis includes 16 samples from three anonymized donors.

The workflow covers Cell Ranger processing, quality control, SCTransform normalization, Harmony integration, clustering, marker identification, cell-type annotation, visualization and enrichment analysis.



## Repository Structure

```text
Human-GW22-24-Cortex-Oligodendroglial-scRNAseq/
├── Bash/
│   └── 00_CellRanger/
│       ├── Donor01_All_CellRanger.sh
│       ├── Donor01_PDGFRA_CellRanger.sh
│       ├── Donor01_EGFR_1_CellRanger.sh
│       ├── Donor01_EGFR_2_CellRanger.sh
│       ├── Donor01_LIFR_1_CellRanger.sh
│       ├── Donor01_LIFR_2_CellRanger.sh
│       ├── Donor02_All_CellRanger.sh
│       ├── Donor02_PDGFRA_CellRanger.sh
│       ├── Donor02_EGFR_CellRanger.sh
│       ├── Donor02_LIFR_1_CellRanger.sh
│       ├── Donor02_LIFR_2_CellRanger.sh
│       ├── Donor02_LIFR_3_CellRanger.sh
│       ├── Donor03_PDGFRA_CellRanger.sh
│       ├── Donor03_EGFR_CellRanger.sh
│       ├── Donor03_LIFR_1_CellRanger.sh
│       └── Donor03_LIFR_2_CellRanger.sh
├── R/
│   ├── 01_QC.R
│   ├── 02_Integration.R
│   ├── 03_Markers.R
│   ├── 04_Annotation.R
│   ├── 05_Figures.R
│   └── 06_Enrichment.R
├── LICENSE
└── README.md
```



Scripts in the `Bash` and `R` directories are numbered according to the recommended analysis workflow.

The Cell Ranger scripts process the 16 individual samples. The R scripts perform quality control, data integration, clustering, marker analysis, cell-type annotation, figure generation and enrichment comparison.



## Data availability

The raw sequencing data are available from Science Data Bank:

[A Gliogenic-lineage Enriched Single-cell Transcriptomic Dataset of Developing Human Cortex](https://doi.org/10.57760/sciencedb.017fa)

Large Cell Ranger outputs and intermediate Seurat objects are not included in this code repository.



## Reproducibility

The scripts in this repository reflect the analysis workflow used for the associated dataset.

Paths to FASTQ files, Cell Ranger references, large input datasets and computing environments must be adjusted before execution on another system. Parameters affecting the reported analyses are retained in the corresponding scripts.



## License

The code in this repository is released under the MIT License.

The sequencing data are distributed separately through Science Data Bank and are not covered by the software license of this repository.
