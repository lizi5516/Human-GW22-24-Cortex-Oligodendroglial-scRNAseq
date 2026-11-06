#!/bin/bash
#SBATCH --job-name=Donor03_EGFR_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor03_EGFR" \
  --fastqs="${FASTQ_ROOT}/Donor03_EGFR" \
  --sample="Donor03_EGFR" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

