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

OUT=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen
VCF=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/05_MAF/joint/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.MAF005.recode.newnames.newchr.vcf
SAMPLES=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/05_MAF/by_pop
#################
# calculate Fst
#################

# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --out $OUT/Fst/LAOvsSTLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --fst-window-size 100000 --fst-window-step 50000 --out $OUT/Fst/LAOvsSTLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep

# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_LAEsamples.txt --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --out $OUT/Fst/LAEvsLAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_LAEsamples.txt --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --fst-window-size 100000 --fst-window-step 50000 --out $OUT/Fst/LAEvsLAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep

# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLO3samples.txt --out $OUT/Fst/STLOvsSTLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLO3samples.txt --fst-window-size 100000 --fst-window-step 50000 --out $OUT/Fst/STLOvsSTLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep

# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLEsamples.txt --out $OUT/Fst/STLOvsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLEsamples.txt --fst-window-size 100000 --fst-window-step 50000 --out $OUT/Fst/STLOvsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep

# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_STLO3samples.txt --weir-fst-pop $SAMPLES/filtered_STLEsamples.txt --out $OUT/Fst/STLO3vsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_STLO3samples.txt --weir-fst-pop $SAMPLES/filtered_STLEsamples.txt --fst-window-size 100000 --fst-window-step 50000 --out $OUT/Fst/STLO3vsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep

# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLEsamples.txt --out $OUT/Fst/LAOvsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLEsamples.txt --fst-window-size 100000 --fst-window-step 50000 --out $OUT/Fst/LAOvsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep

# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --out $OUT/Fst/LAOvsSTLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed
# vcftools --vcf $VCF --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --fst-window-size 100000 --fst-window-step 50000 --out $OUT/Fst/LAOvsSTLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin50KBstep


vcftools --vcf $VCF --chr 10 --weir-fst-pop $SAMPLES/filtered_LAEsamples.txt --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --fst-window-size 100000 --fst-window-step 100000 --out $OUT/Fst/LAEvsLAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin.nooverlap

vcftools --vcf $VCF --chr 10 --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLO3samples.txt --fst-window-size 100000 --fst-window-step 100000 --out $OUT/Fst/STLOvsSTLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin.nooverlap

vcftools --vcf $VCF --chr 10 --weir-fst-pop $SAMPLES/filtered_STLOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLEsamples.txt --fst-window-size 100000 --fst-window-step 100000 --out $OUT/Fst/STLOvsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin.nooverlap

vcftools --vcf $VCF --chr 10 --weir-fst-pop $SAMPLES/filtered_STLO3samples.txt --weir-fst-pop $SAMPLES/filtered_STLEsamples.txt --fst-window-size 100000 --fst-window-step 100000 --out $OUT/Fst/STLO3vsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin.nooverlap

vcftools --vcf $VCF --chr 10 --weir-fst-pop $SAMPLES/filtered_LAOsamples.txt --weir-fst-pop $SAMPLES/filtered_STLEsamples.txt --fst-window-size 100000 --fst-window-step 100000 --out $OUT/Fst/LAOvsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.100KBwin.nooverlap




#################
# calculate TajD
#################
# vcftools --vcf $VCF --TajimaD 10000 --keep $SAMPLES/filtered_STLOsamples.txt --out $OUT/TajD/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf $VCF --TajimaD 10000 --keep $SAMPLES/filtered_STLO3samples.txt --out $OUT/TajD/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf $VCF --TajimaD 10000 --keep $SAMPLES/filtered_STLEsamples.txt --out $OUT/TajD/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf $VCF --TajimaD 10000 --keep $SAMPLES/filtered_LAOsamples.txt --out $OUT/TajD/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf $VCF --TajimaD 10000 --keep $SAMPLES/filtered_LAEsamples.txt --out $OUT/TajD/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB

# vcftools --vcf $VCF --TajimaD 100000 --keep $SAMPLES/filtered_STLOsamples.txt --out $OUT/TajD/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD100KB

# vcftools --vcf $VCF --TajimaD 100000 --keep $SAMPLES/filtered_STLO3samples.txt --out $OUT/TajD/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD100KB

# vcftools --vcf $VCF --TajimaD 100000 --keep $SAMPLES/filtered_STLEsamples.txt --out $OUT/TajD/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD100KB

# vcftools --vcf $VCF --TajimaD 100000 --keep $SAMPLES/filtered_LAOsamples.txt --out $OUT/TajD/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD100KB

# vcftools --vcf $VCF --TajimaD 100000 --keep $SAMPLES/filtered_LAEsamples.txt --out $OUT/TajD/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD100KB


# vcftools --vcf $VCF --TajimaD 1000000 --keep $SAMPLES/filtered_STLOsamples.txt --out $OUT/TajD/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD1MB

# vcftools --vcf $VCF --TajimaD 1000000 --keep $SAMPLES/filtered_STLO3samples.txt --out $OUT/TajD/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD1MB

# vcftools --vcf $VCF --TajimaD 1000000 --keep $SAMPLES/filtered_STLEsamples.txt --out $OUT/TajD/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD1MB

# vcftools --vcf $VCF --TajimaD 1000000 --keep $SAMPLES/filtered_LAOsamples.txt --out $OUT/TajD/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD1MB

# vcftools --vcf $VCF --TajimaD 1000000 --keep $SAMPLES/filtered_LAEsamples.txt --out $OUT/TajD/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD1MB

# vcftools --vcf $VCF --TajimaD 5000000 --keep $SAMPLES/filtered_STLOsamples.txt --out $OUT/TajD/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD5MB

# vcftools --vcf $VCF --TajimaD 5000000 --keep $SAMPLES/filtered_STLO3samples.txt --out $OUT/TajD/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD5MB

# vcftools --vcf $VCF --TajimaD 5000000 --keep $SAMPLES/filtered_STLEsamples.txt --out $OUT/TajD/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD5MB

# vcftools --vcf $VCF --TajimaD 5000000 --keep $SAMPLES/filtered_LAOsamples.txt --out $OUT/TajD/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD5MB

# vcftools --vcf $VCF --TajimaD 5000000 --keep $SAMPLES/filtered_LAEsamples.txt --out $OUT/TajD/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD5MB

#################
# calculate Pi
#################
# vcftools --vcf $VCF --window-pi 10000 --keep $SAMPLES/filtered_STLOsamples.txt --out $OUT/pi/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi10KB

# vcftools --vcf $VCF --window-pi 10000 --keep $SAMPLES/filtered_STLO3samples.txt --out $OUT/pi/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi10KB

# vcftools --vcf $VCF --window-pi 10000 --keep $SAMPLES/filtered_STLEsamples.txt --out $OUT/pi/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi10KB

# vcftools --vcf $VCF --window-pi 10000 --keep $SAMPLES/filtered_LAOsamples.txt --out $OUT/pi/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi10KB

# vcftools --vcf $VCF --window-pi 10000 --keep $SAMPLES/filtered_LAEsamples.txt --out $OUT/pi/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi10KB

# vcftools --vcf $VCF --window-pi 100000 --keep $SAMPLES/filtered_STLOsamples.txt --out $OUT/pi/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi100KB

# vcftools --vcf $VCF --window-pi 100000 --keep $SAMPLES/filtered_STLO3samples.txt --out $OUT/pi/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi100KB

# vcftools --vcf $VCF --window-pi 100000 --keep $SAMPLES/filtered_STLEsamples.txt --out $OUT/pi/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi100KB

# vcftools --vcf $VCF --window-pi 100000 --keep $SAMPLES/filtered_LAOsamples.txt --out $OUT/pi/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi100KB

# vcftools --vcf $VCF --window-pi 100000 --keep $SAMPLES/filtered_LAEsamples.txt --out $OUT/pi/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi100KB

# vcftools --vcf $VCF --window-pi 1000000 --keep $SAMPLES/filtered_STLOsamples.txt --out $OUT/pi/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi1MB

# vcftools --vcf $VCF --window-pi 1000000 --keep $SAMPLES/filtered_STLO3samples.txt --out $OUT/pi/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi1MB

# vcftools --vcf $VCF --window-pi 1000000 --keep $SAMPLES/filtered_STLEsamples.txt --out $OUT/pi/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi1MB

# vcftools --vcf $VCF --window-pi 1000000 --keep $SAMPLES/filtered_LAOsamples.txt --out $OUT/pi/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi1MB

# vcftools --vcf $VCF --window-pi 1000000 --keep $SAMPLES/filtered_LAEsamples.txt --out $OUT/pi/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi1MB

# vcftools --vcf $VCF --window-pi 5000000 --keep $SAMPLES/filtered_STLOsamples.txt --out $OUT/pi/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB

# vcftools --vcf $VCF --window-pi 5000000 --keep $SAMPLES/filtered_STLO3samples.txt --out $OUT/pi/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB

# vcftools --vcf $VCF --window-pi 5000000 --keep $SAMPLES/filtered_STLEsamples.txt --out $OUT/pi/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB

# vcftools --vcf $VCF --window-pi 5000000 --keep $SAMPLES/filtered_LAOsamples.txt --out $OUT/pi/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB

# vcftools --vcf $VCF --window-pi 5000000 --keep $SAMPLES/filtered_LAEsamples.txt --out $OUT/pi/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB

#################
# calculate relatedness and het
#################

# vcftools --vcf GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.vcf --relatedness2 --out GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed

# vcftools --vcf $VCF --het --out $OUT/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed

#################
# create new vcf with just outliers (and samples for LAOvSTLO)
#################

# vcftools --vcf $VCF --keep $SAMPLES/LAO_STLO_samples.txt --bed $OUT/Fst/outliers/LAOvLAEoutlier_windows.bed --out $OUT/Fst/outliers/LAOvLAEoutlier_window --recode --keep-INFO-all

# vctools --vcf ./LAOvLAEoutlier_window.recode.vcf -keep ../../../../09_filterVCFs/05_MAF/by_pop/filtered_STLOsamples.txt --out

#################
# calculate LD chr10
#################

# vcftools --vcf $VCF --geno-r2 --chr 10 --thin 50000 --keep $SAMPLES/filtered_STLOsamples.txt --out $OUT/linkage/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.thin50K

# vcftools --vcf $VCF --geno-r2 --chr 10 --thin 50000 --keep $SAMPLES/filtered_STLO3samples.txt --out $OUT/linkage/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.thin50K

# vcftools --vcf $VCF --geno-r2 --chr 10 --thin 50000 --keep $SAMPLES/filtered_STLEsamples.txt --out $OUT/linkage/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.thin50K

# vcftools --vcf $VCF --geno-r2 --chr 10 --thin 50000 --keep $SAMPLES/filtered_LAOsamples.txt --out $OUT/linkage/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.thin50K

# vcftools --vcf $VCF --geno-r2 --chr 10 --thin 50000 --keep $SAMPLES/filtered_LAEsamples.txt --out $OUT/linkage/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.thin50K

#################
# look at site depth for chr10
#################

# vcftools --vcf $VCF --chr 10 --keep $SAMPLES/filtered_STLOsamples.txt --site-mean-depth --out $OUT/STLO_chr10

# vcftools --vcf $VCF --chr 10 --keep $SAMPLES/filtered_STLO3samples.txt --site-mean-depth --out $OUT/STLO3_chr10

# vcftools --vcf $VCF --chr 10 --keep $SAMPLES/filtered_STLEsamples.txt --site-mean-depth --out $OUT/STLE_chr10

# vcftools --vcf $VCF --chr 10 --keep $SAMPLES/filtered_LAOsamples.txt --site-mean-depth --out $OUT/LAO_chr10

# vcftools --vcf $VCF --chr 10 --keep $SAMPLES/filtered_LAEsamples.txt --site-mean-depth --out $OUT/LAE_chr10

#################
# make a vcf with all sampels just for chr10
#################

# vcftools --vcf $VCF --chr 10 --out $OUT/GL_Pink_Salmon_chr10 --recode --keep-INFO-all
