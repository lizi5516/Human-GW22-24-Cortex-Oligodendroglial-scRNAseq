#!/bin/bash
#SBATCH --job-name=Donor02_LIFR_1_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor02_LIFR_1" \
  --fastqs="${FASTQ_ROOT}/Donor02_LIFR_1" \
  --sample="Donor02_LIFR_1" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

