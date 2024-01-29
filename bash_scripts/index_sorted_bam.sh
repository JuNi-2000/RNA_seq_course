#!/bin/bash
#SBATCH --array=1-16
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=feature_counts
#SBATCH --output=bam_sort_%J.out
#SBATCH --error=bam_sort_%J.err
#SBATCH --partition=pall


# define variables
WORKDIR="/data/users/jniklaus2/RNA_seq"
OUTDIR="$WORKDIR/mapping"
SAMPLELIST="$WORKDIR/metadata/sample_list_bam_index.tsv"

SAMPLENAME=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
SAMPLEPATH=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`

cd $WORKDIR
module load UHTS/Analysis/samtools/0.1.19
samtools index $SAMPLEPATH