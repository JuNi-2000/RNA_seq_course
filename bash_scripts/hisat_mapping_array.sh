#!/bin/bash
#SBATCH --array=1-16
#SBATCH --time=03:00:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=4
#SBATCH --job-name=HISAT2_RNA_Mapping
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pall

# define variables
WORKDIR="/data/users/jniklaus2/RNA_seq"
OUTDIR="$WORKDIR/results"
SAMPLELIST="$WORKDIR/metadata/samplelist.tsv"

SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

OUTFILE="$OUTDIR/${SAMPLE}.sam"

############################

mkdir -p $OUTDIR
module load UHTS/Aligner/hisat/2.2.1
hisat2 -q --rna-strandness RF -x $WORKDIR/index_hisat/genome -1 $READ1 -2 $READ2 -S $OUTFILE

#Testing wheter array finds the right files
echo "Run task for $SAMPLE with $READ1 and $READ2; corresponding Task ID: $SLURM_ARRAY_TASK_ID" > $OUTFILE
