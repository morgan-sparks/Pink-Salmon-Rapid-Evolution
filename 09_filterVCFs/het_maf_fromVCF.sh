#!/bin/bash

#SBATCH --job-name=het_maf_fromVCF
#SBATCH -A beagle
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sparks35@purdue.edu

module purge

module load bioinfo
module load r

FILEPATH=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/09_filterVCFs/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/gntyp_missingness/MAF_filter/MAF_0.2/

Rscript $FILEPATH/heterozygosity_maf_from_vcf_FAST_v2.R

FILEPATH=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/09_filterVCFs/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/gntyp_missingness/MAF_filter/MAF_0.3/

Rscript $FILEPATH/heterozygosity_maf_from_vcf_FAST_v2.R