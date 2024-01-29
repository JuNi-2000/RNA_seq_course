#!/bin/bash
#SBATCH --array=1-1
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4
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


OUTFILE="$OUTDIR/${SAMPLENAME}_counts"
#cho "Run task for $SAMPLENAME with $SAMPLEPATH" > $OUTFILE

cd $WORKDIR
module load UHTS/Analysis/subread/2.0.1
featureCounts -S 2 -a /data/users/jniklaus2/RNA_seq/reference_genome/Mus_musculus.GRCm39.110.gtf -o $WORKDIR/feature_counts/counts.txt \
/data/users/jniklaus2/RNA_seq/mapping/sorted_bam/1918_sorted.bam

#-S2 is the strandedness, -a is the location of annotation file





module load UHTS/Analysis/samtools/0.1.19
samtools sort $SAMPLEPATH $OUTFILE
#do not use -o for output! That was used in older versions, now it leads to stout printing