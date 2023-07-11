#!/bin/bash

#SBATCH --job-name=align_LAO_199_run2
#SBATCH -A standby
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 126
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge
module load bioinfo
module load bwa/0.7.17

PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon

cd $PROJHOME/data/seqs/aligned_reads_Ogor1.0/

bwa mem -t 126 -M $PROJHOME/data/assemblies/Ogor_1.0/GCA_017355495.1_Ogor_1.0_genomic.fna $PROJHOME/data/seqs/trimmed_seqs/trimmed_fastq_v2/trimmed_paired_LAO_199_R1_run2.fastq $PROJHOME/data/seqs/trimmed_seqs/trimmed_fastq_v2/trimmed_paired_LAO_199_R2_run2.fastq > $PROJHOME/data/seqs/aligned_reads_Ogor1.0/02_batch_output/alignment_v2/LAO_199_run2_Ogor1.0_aligned.bam
