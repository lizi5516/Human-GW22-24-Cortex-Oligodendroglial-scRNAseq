#!/bin/bash
#SBATCH --job-name=Donor02_LIFR_3_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor02_LIFR_3" \
  --fastqs="${FASTQ_ROOT}/Donor02_LIFR_3" \
  --sample="Donor02_LIFR_3" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

