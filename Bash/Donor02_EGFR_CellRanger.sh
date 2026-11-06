#!/bin/bash
#SBATCH --job-name=Donor02_EGFR_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor02_EGFR" \
  --fastqs="${FASTQ_ROOT}/Donor02_EGFR" \
  --sample="Donor02_EGFR" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

