#!/bin/bash
#SBATCH --job-name=Donor01_PDGFRA_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor01_PDGFRA" \
  --fastqs="${FASTQ_ROOT}/Donor01_PDGFRA" \
  --sample="Donor01_PDGFRA" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

