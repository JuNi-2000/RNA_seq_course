#!/bin/bash
#SBATCH --array=1-16
#SBATCH --time=01:00:00
#SBATCH --mem=25g
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sort_bam
#SBATCH --output=bam_sort_%J.out
#SBATCH --error=bam_sort_%J.err
#SBATCH --partition=pall

# define variables
WORKDIR="/data/users/jniklaus2/RNA_seq"
OUTDIR="$WORKDIR/mapping"
SAMPLELIST="$WORKDIR/metadata/sample_names_bamfiles.tsv"

SAMPLENAME=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
SAMPLEPATH=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`


OUTFILE="$OUTDIR/${SAMPLENAME}_sorted.bam"
#cho "Run task for $SAMPLENAME with $SAMPLEPATH" > $OUTFILE

cd $WORKDIR
module load UHTS/Analysis/samtools/0.1.19
samtools sort $SAMPLEPATH $OUTFILE
#do not use -o for output! That was used in older versions, now it leads to stout printing