#!/bin/bash
#SBATCH --job-name=filter_invariant_sites
#SBATCH -A beagle
#SBATCH -t 300:00:00
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge
module load bioinfo
module load vcftools/0.1.14

OUT=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/04_gntypMissing/invariant_filter
VCF=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/04_gntypMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.vcf


# filtering by pop and maf 0.03333333 which is 1 variant/30 samples (i.e., get rid of all invariant sites)

vcftools --vcf $VCF --keep $OUT/filtered_STLOsamples.txt --maf 0.03333333 --recode --recode-INFO-all --out $OUT/STLO_noMAF.filt

vcftools --vcf $VCF --keep $OUT/filtered_STLO3samples.txt --maf 0.03333333 --recode --recode-INFO-all --out $OUT/STLO3_noMAF.filt

vcftools --vcf $VCF --keep $OUT/filtered_STLEsamples.txt --maf 0.03333333 --recode --recode-INFO-all --out $OUT/STLE_noMAF.filt

vcftools --vcf $VCF --keep $OUT/filtered_LAOsamples.txt --maf 0.03333333 --recode --recode-INFO-all --out $OUT/LAO_noMAF.filt

vcftools --vcf $VCF --keep $OUT/filtered_LAEsamples.txt --maf 0.03333333 --recode --recode-INFO-all --out $OUT/LAE_noMAF.filt