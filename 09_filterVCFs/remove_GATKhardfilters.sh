#!/bin/bash

#SBATCH --job-name=remove_GATK_hardfilters
#SBATCH -A beagle
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sparks35@purdue.edu

module purge

module load bioinfo
module load vcftools/0.1.16

FILTERING=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/01_GATK_hardfilter 

vcftools --gzvcf $FILTERING/GL_PinkSalmon_even_GATKfilter_snps.vcf.gz --remove-filtered-all --recode --recode-INFO-all --out \
$FILTERING/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.vcf.gz

