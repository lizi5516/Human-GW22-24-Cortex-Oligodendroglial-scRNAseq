#!/bin/bash
#SBATCH --job-name=Donor03_PDGFRA_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor03_PDGFRA" \
  --fastqs="${FASTQ_ROOT}/Donor03_PDGFRA" \
  --sample="Donor03_PDGFRA" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

