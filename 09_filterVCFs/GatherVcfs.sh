#!/bin/bash

#SBATCH --job-name=gather_chroms
#SBATCH -A beagle
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge
PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
FILTERED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/01_GATK_hardfilter


## note: I leave out NC_060199.1_even_GATKfilter_snps.vcf.gz because that is the Y-chrom scafffold
# Below commnad uses GatherVcfs to merge chromosomes into a single vcf file with all jointly called
# individuals. These are GATK marked (not yet filtered) SNPs only.

/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx40g " GatherVcfs \
-I $FILTERED/NC_060173.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060174.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060175.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060176.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060177.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060178.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060179.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060180.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060181.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060182.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060183.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060184.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060185.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060186.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060187.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060188.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060189.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060190.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060191.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060192.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060193.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060194.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060195.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060196.1_even_GATKfilter_snps.vcf.gz \
-I $FILTERED/NC_060197.1_even_GATKfilter_snps.vcf.gz \
-O $FILTERED/GL_PinkSalmon_even_GATKfilter_allsnps.vcf.gz



