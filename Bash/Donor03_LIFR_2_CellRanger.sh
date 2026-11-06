#!/bin/bash
#SBATCH --job-name=Donor03_LIFR_2_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor03_LIFR_2" \
  --fastqs="${FASTQ_ROOT}/Donor03_LIFR_2" \
  --sample="Donor03_LIFR_2" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

