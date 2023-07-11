#!/bin/bash
#SBATCH --job-name=calc_PopGen
#SBATCH -A beagle
#SBATCH -t 300:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge 
module load bioinfo
module load vcftools/0.1.14

cd /scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/09_filterVCFs/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/gntyp_missingness/MAF_filter/MAF_0.05

#################
# calculate Fst
#################

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --weir-fst-pop ../MAF_byPop/filtered_LAOsamples.txt --weir-fst-pop ../MAF_byPop/filtered_STLOsamples.txt --out LAOvsSTLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --weir-fst-pop ../MAF_byPop/filtered_LAOsamples.txt --weir-fst-pop ../MAF_byPop/filtered_STLOsamples.txt --fst-window-size 100000 --fst-window-step 50000 --out LAOvsSTLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.vcf --weir-fst-pop ../MAF_byPop/filtered_LAEsamples.txt --weir-fst-pop ../MAF_byPop/filtered_LAOsamples.txt --out LAEvsLAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.vcf --weir-fst-pop ../MAF_byPop/filtered_LAEsamples.txt --weir-fst-pop ../MAF_byPop/filtered_LAOsamples.txt --fst-window-size 100000 --fst-window-step 50000 --out LAEvsLAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --weir-fst-pop ../MAF_byPop/filtered_STLOsamples.txt --weir-fst-pop ../MAF_byPop/filtered_STLO3samples.txt --out STLOvsSTLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --weir-fst-pop ../MAF_byPop/filtered_STLOsamples.txt --weir-fst-pop ../MAF_byPop/filtered_STLO3samples.txt --fst-window-size 100000 --fst-window-step 50000 --out STLOvsSTLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep

#################
# calculate TajD
#################
# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --TajimaD 10000 --keep ../MAF_byPop/filtered_STLOsamples.txt --out STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --TajimaD 10000 --keep ../MAF_byPop/filtered_STLO3samples.txt --out STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --TajimaD 10000 --keep ../MAF_byPop/filtered_STLEsamples.txt --out STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --TajimaD 10000 --keep ../MAF_byPop/filtered_LAOsamples.txt --out LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --TajimaD 100000 --keep ../MAF_byPop/filtered_STLOsamples.txt --out STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --TajimaD 100000 --keep ../MAF_byPop/filtered_STLO3samples.txt --out STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --TajimaD 100000 --keep ../MAF_byPop/filtered_STLEsamples.txt --out STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf --TajimaD 100000 --keep ../MAF_byPop/filtered_LAOsamples.txt --out LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

#################
# calculate relatedness and het
#################

vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.vcf --relatedness2 --out GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed

vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.vcf --relatedness2 --het GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
