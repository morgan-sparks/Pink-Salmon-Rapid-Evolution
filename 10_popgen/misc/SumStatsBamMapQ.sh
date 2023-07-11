#!/bin/bash
#SBATCH -n 15
#SBATCH -A beagle
#SBATCH -t 300:00:00
#SBATCH --mail-user=spark35@purdue.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=SumStatsBamMapQ


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
 samtools view $file | awk '{sum+=$5} END { print sum/NR}' >> $SCRIPTS/BamAvgMapQ.txt 
 samtools view $file | awk '{sum+=$5; sumsq+=$5*$5} END {print sqrt(sumsq/NR - (sum/NR)**2)}' >> $SCRIPTS/BamStdevMapq.txt ; done
 #echo $file >> $SCRIPTS/combinedBamsFileName.txt ; done ###this file was already created in a different process


echo "END TIME: `date`"

 
