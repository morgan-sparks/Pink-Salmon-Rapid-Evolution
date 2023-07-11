#!/bin/bash

#SBATCH --job-name=hardfilter_CM029860.1
#SBATCH -A beagle
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge

PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
CALLSNPS=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/08_genotypeGVCFs/run2
FILTERED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/09_filterVCFs/run2_allsamps/gatk_HardFilter


#select SNP variants only
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx20g " SelectVariants \
-V $CALLSNPS/allVariants/CM029860.1.vcf.gz \
-select-type SNP \
-O $CALLSNPS/SNPsOnly/CM029860.1_snps.vcf.gz 

#filter using GATK recommended hard filter thresholds https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx20g " VariantFiltration \
-V $CALLSNPS/SNPsOnly/CM029860.1_snps.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O $FILTERED/CM029860.1_hardfiltered_snps.vcf.gz
