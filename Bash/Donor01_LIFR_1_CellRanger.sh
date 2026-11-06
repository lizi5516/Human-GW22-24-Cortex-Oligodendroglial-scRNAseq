#!/bin/bash
#SBATCH --job-name=Donor01_LIFR_1_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor01_LIFR_1" \
  --fastqs="${FASTQ_ROOT}/Donor01_LIFR_1" \
  --sample="Donor01_LIFR_1" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

