#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800M
#SBATCH --time=00:10:00
#SBATCH --job-name=multiqc

module load UHTS/Analysis/MultiQC/1.8
cd /data/users/jniklaus2/RNA_seq/fastqc
multiqc .