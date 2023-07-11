#!/bin/bash

#SBATCH --job-name=HWEfromVCF
#SBATCH -A beagle
#SBATCH -t 330:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sparks35@purdue.edu

module purge

module load bioinfo
module load r

cd /scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/09_filterVCFs/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/gntyp_missingness/MAF_filter/MAF_0.05

FILEPATH=/scratch/bell/sparks35/GL_Pink_Salmon/scripts/09_filterVCFs

Rscript $FILEPATH/HWE_from_VCF.R

