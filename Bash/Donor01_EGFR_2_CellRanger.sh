#!/bin/bash
#SBATCH --job-name=Donor01_EGFR_2_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor01_EGFR_2" \
  --fastqs="${FASTQ_ROOT}/Donor01_EGFR_2" \
  --sample="Donor01_EGFR_2" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

