#!/bin/bash
#SBATCH --job-name=Donor02_All_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor02_All" \
  --fastqs="${FASTQ_ROOT}/Donor02_All" \
  --sample="Donor02_All" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

