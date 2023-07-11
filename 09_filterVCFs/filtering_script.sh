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

# vcftools --gzvcf ./01_GATK_hardfilter/GL_PinkSalmon_even_GATKfilter_snps.vcf.gz --remove-filtered-all \
# --recode --recode-INFO-all --out .//01_GATK_hardfilter/GL_PinkSalmon_even_GATKfilter.GATKrem.snps

#  vcftools --vcf ./01_GATK_hardfilter/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.vcf --min-alleles 2 --max-alleles 2 \
#  --recode --recode-INFO "biallelic" --out ./02_remove_nonbiallelic/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.

## script for findng indv w/ >20% missingness

# vcftools --vcf ./02_remove_nonbiallelic/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.vcf --missing-indv \
# --out ./03_indvMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.vcf

## script to remove those individuals

# vcftools --vcf ./02_remove_nonbiallelic/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.vcf \
# --remove ./03_indvMissing/indv80.txt \
# --recode --recode-INFO-all --out ./03_indvMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80

## script to remove genotypes w/ > 20% missingness

# vcftools --vcf ./03_indvMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.vcf \
# --max-missing-count 52 --recode --recode-INFO-all --out ./04_gntypMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80

## filter jointly on MAF 0.05
s
# vcftools --vcf ./04_gntypMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.vcf \
# --maf .05 --recode --recode-INFO-all \
# --out ./05_MAF/joint/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.MAF005

## Filter by pop on MAF 0.05

for FILE in LAE LAO STLE STLO3 STLO
do
vcftools --vcf ./04_gntypMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.vcf \
--keep ./05_MAF/by_pop/filtered_${FILE}samples.txt \
--maf .05 --recode --recode-INFO-all \
--out ./05_MAF/by_pop/${FILE}samples_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.MAF005
done


