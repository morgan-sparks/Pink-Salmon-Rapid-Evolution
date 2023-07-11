#!/bin/bash
#SBATCH -n 6
#SBATCH -A beagle
#SBATCH -t 30:00:00
#SBATCH --mail-user=spark35@purdue.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=BamSummary


	###Single###
	###Alignments ###
#####################################
### Programs Used ###################
#####################################
module load bioinfo
module load samtools

echo "START TIME: `date`"

DATA=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/04_merged_bams
SCRIPTS=/scratch/bell/sparks35/GL_Pink_Salmon/scripts/10_popgen/misc
cd $DATA/

echo ""

echo ""

for file in *bam; do
 samtools view -c -q 20 $file >> $SCRIPTS/combinedBamsQ20Reads.txt 
 samtools view -c -q 60 $file >> $SCRIPTS/combinedBamsQ60Reads.txt ; done


echo "END TIME: `date`"

 
