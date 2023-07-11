#!/bin/bash

#SBATCH --job-name=addreadgroups_LAO_195_run2
#SBATCH -A standby
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge
module load bioinfo
module load picard-tools/2.20.8

PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
ALIGNED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/02_aligned_reads/alignment_even
RGADD=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/03_addedRG_bam/alignment_even


java -Xmx40g -jar /group/bioinfo/apps/apps/picard-tools-2.9.0/picard.jar AddOrReplaceReadGroups \
I=$ALIGNED/LAO_195_run2_Ogor1.0_aligned.bam \
O=$RGADD/LAO_195_run2_Ogor1.0_addedRG.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=LAO_195 RGPU=run2
