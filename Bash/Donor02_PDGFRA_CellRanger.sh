#!/bin/bash
#SBATCH --job-name=Donor02_PDGFRA_CR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

cellranger count \
  --id="Donor02_PDGFRA" \
  --fastqs="${FASTQ_ROOT}/Donor02_PDGFRA" \
  --sample="Donor02_PDGFRA" \
  --transcriptome="${CELLRANGER_TRANSCRIPTOME}" \
  --create-bam=true

