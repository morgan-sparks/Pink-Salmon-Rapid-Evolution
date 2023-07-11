#!/bin/bash

#SBATCH --job-name=mergebams_STL_033
#SBATCH -A standby
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge
module load bioinfo
module load picard-tools/2.20.8
module load samtools/1.12

PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
ADDRG=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0//03_addedRG_bam/alignment_even
MERGED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/04_merged_bams/alignment_even


java -Xmx124g -jar /group/bioinfo/apps/apps/picard-tools-2.9.0/picard.jar MergeSamFiles \
I=$ADDRG/STL_033_run1_Ogor1.0_addedRG.bam \
I=$ADDRG/STL_033_run2_Ogor1.0_addedRG.bam \
SORT_ORDER=coordinate \
O=$MERGED/STL_033_Ogor1.0_merged.bam

samtools flagstat -@ 126 -O tsv $MERGED/STL_033_Ogor1.0_merged.bam > $MERGED/merged_flagstat/STL_033_rgroups_flagstat_out.txt
