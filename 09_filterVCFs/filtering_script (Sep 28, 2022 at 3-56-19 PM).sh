#!/bin/bash

#SBATCH --job-name=filtering
#SBATCH -A beagle
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sparks35@purdue.edu

module purge

module load bioinfo
module load vcftools/0.1.16

cd /scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs

 vcftools --vcf ./01_GATK_hardfilter/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.vcf --min-alleles 2 --max-alleles 2 \
 --recode --recode-INFO "biallelic" --out ./02_remove_nonbiallelic/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.