#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=01:30:00
#SBATCH --mem=4g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=sam_to_bam_conversion
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pall

# define variables
WORKDIR="/data/users/jniklaus2/RNA_seq"
OUTDIR="$WORKDIR/mapping"
SAMPLELIST="$WORKDIR/metadata/sample_names_samfiles.tsv"

SAMPLENAME=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
SAMPLEPATH=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`


OUTFILE="$OUTDIR/${SAMPLENAME}.txt"
#cho "Run task for $SAMPLENAME with $SAMPLEPATH" > $OUTFILE

cd WORKDIR
module load UHTS/Analysis/samtools/0.1.19
samtools view -bS $SAMPLEPATH > $SAMPLENAME.bam