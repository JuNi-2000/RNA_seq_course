#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800M
#SBATCH --time=01:00:00
#SBATCH --job-name=hisat_index

cd /data/users/jniklaus2/RNA_seq
module load UHTS/Aligner/hisat/2.2.1

hisat2-build -p 16 Mus_musculus.GRCm39.dna.primary_assembly.fa genome


