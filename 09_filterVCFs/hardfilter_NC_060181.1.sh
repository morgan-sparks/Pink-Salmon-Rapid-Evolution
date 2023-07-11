#!/bin/bash

#SBATCH --job-name=hardfilter_NC_060181.1
#SBATCH -A beagle
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge

PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
CALLSNPS=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/08_genotypeGVCFs
FILTERED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/01_GATK_hardfilter


#select SNP variants only
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx20g " SelectVariants \
-V $CALLSNPS/allVariants/NC_060181.1.even.vcf.gz \
-select-type SNP \
-O $CALLSNPS/SNPsOnly/NC_060181.1_snps.vcf.gz 

#filter using GATK recommended hard filter thresholds https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx20g " VariantFiltration \
-V $CALLSNPS/SNPsOnly/NC_060181.1_snps.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O $FILTERED/NC_060181.1_even_GATKfilter_snps.vcf.gz
