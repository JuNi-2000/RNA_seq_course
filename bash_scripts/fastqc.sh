#!/bin/bash
#SBATCH --array=1-32
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=FAST_qc
#SBATCH --error=array_%J.err
#SBATCH --partition=pall


WORKDIR="/data/users/jniklaus2/RNA_seq/"

SAMPLELIST="$WORKDIR/metadata/sample_list_fastqc.tsv"

SAMPLENAME=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
SAMPLEPATH=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`


cd $WORKDIR
module add UHTS/Quality_control/fastqc/0.11.9
fastqc --extract /data/users/jniklaus2/RNA_seq/reads/reads/$SAMPLENAME -o /data/users/jniklaus2/RNA_seq/fastqc

