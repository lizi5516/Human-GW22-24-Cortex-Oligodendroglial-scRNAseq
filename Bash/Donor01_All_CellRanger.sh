#!/bin/bash

#SBATCH --job-name=Donor01_All_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32


cellranger count \
  --id="Donor01_All" \
  --fastqs="${FASTQ_ROOT}/Donor01_All" \
  --sample="Donor01_All" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true
