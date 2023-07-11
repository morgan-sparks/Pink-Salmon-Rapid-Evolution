#!/bin/bash
#SBATCH --job-name=snpEff_snpSift
#SBATCH -A mcclintock
#SBATCH -t 300:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge
module load biocontainers
module load snpeff

cd /scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers

### annotate with snpEff

## Zfst sites
snpEff ann -c ~/snpEff/snpEff/snpEff.config OgorEven_v1.0 -Xms25g STLO_Zfst_windows.recode.vcf > ./snpEff_results/STLO_results/STLO_Zfst_windows.ann.vcf

snpEff ann -c ~/snpEff/snpEff/snpEff.config OgorEven_v1.0 -Xms25g  LAO_Zfst_windows.recode.vcf > ./snpEff_results/LAO_results/LAO_Zfst_windows.ann.vcf


## inversion sites
snpEff ann -c ~/snpEff/snpEff/snpEff.config OgorEven_v1.0 -Xms25g STLO_chr10_inversion.recode.vcf > ./snpEff_results/STLO_results/STLO_ch10_inversion.ann.vcf

snpEff ann -c ~/snpEff/snpEff/snpEff.config OgorEven_v1.0 -Xms25g LAO_chr10_inversion.recode.vcf > ./snpEff_results/LAO_results/LAO_ch10_inversion.ann.vcf

## ZFst SNPs in ZFst windows
snpEff ann -c ~/snpEff/snpEff/snpEff.config OgorEven_v1.0 -Xms25g STLO_ZFstoutliers_inZFstwindows.recode.vcf > ./snpEff_results/STLO_ZFstoutliers_inZFstwindows.ann.vcf

cat STLO_ZFstoutliers_inZFstwindows.ann.vcf | ~/snpEff/snpEff/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -Xms25g - CHROM POS REF ALT "ANN[*].EFFECT" "EFF[*].GENE" > STLO_snpEff_ZfstSNPs_inZfstwindows.txt

### extract variants and effects

### Zfst sites
cat STLO_Zfst_windows.ann.vcf | ~/snpEff/snpEff/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -Xms25g - CHROM POS REF ALT "ANN[*].EFFECT" "EFF[*].GENE" > STLO_snpEff_Zfst_SNPs.txt

cat LAO_Zfst_windows.ann.vcf | ~/snpEff/snpEff/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -Xms25g - CHROM POS REF ALT "ANN[*].EFFECT" "EFF[*].GENE" > LAO_snpEff_Zfst_SNPs.txt


### inversion sites
SnpSift extractFields -Xms25g STLO_ch10_inversion.ann.vcf CHROM POS REF ALT "ANN[*].EFFECT" "EFF[*].GENE" > STLO_snpEff_SNPs.txt

cat STLO_ch10_inversion.ann.vcf | ~/snpEff/snpEff/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -Xms25g - CHROM POS REF ALT "ANN[*].EFFECT" "EFF[*].GENE" > STLO_snpEff_SNPs.txt

cat LAO_ch10_inversion.ann.vcf | ~/snpEff/snpEff/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -Xms25g - CHROM POS REF ALT "ANN[*].EFFECT" "EFF[*].GENE" > LAO_snpEff_SNPs.txt

## build database
#snpEff build -c ./snpEff.config -gtf22 -v GCF_021184085.1