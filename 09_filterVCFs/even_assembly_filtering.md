# Filtering for Even Assembly samples

These commands were mostly  run with an interactive job in vcftools/1.16.0 or as a generalized bash script in the format below

```
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

### one of the below filtering commands
```
***

## Remove GATK filters

First step is to take chromosome specific vcfs and use GATKs GatherVcfs function to make single VCF called `GL_PinkSalmon_even_GATKfilter_snps.vcf.gz`

You can count the number of SNPs using `grep -v "^#" FILE.vcf | wc -l &`

**Original unfiltered file has 24,312,640 SNPs**

The use `remove_GATKhardfilters.sh` to remove the marked GATK filters from the `hardfilter_NC*.sh` scripts (which run on chromosomes).

```
vcftools --gzvcf ./01_GATK_hardfilter/GL_PinkSalmon_even_GATKfilter_snps.vcf.gz --remove-filtered-all --recode --recode-INFO-all --out ./01_GATK_hardfilter/GL_PinkSalmon_even_GATKfilter.GATKrem.snps
```

**GATK filtering brought file to 20,370,470  SNPs**

***

## Remove non-biallelic sites

Now filter on non-biallelic SNPs

```
 vcftools --vcf ./01_GATK_hardfilter/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.vcf --min-alleles 2 --max-alleles 2  --recode --recode-INFO-all --out ./02_remove_nonbiallelic/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic
```

**Filtering non-bialleic sites brought file to 18,362,895 SNPs**

***

## Find indviduals with with >20% missingness and then remove them

First list individuals with >20% missingness

```
vcftools --vcf ./02_remove_nonbiallelic/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.vcf --missing-indv --out ./03_indvMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.vcf
```

To find which individuals have >20% missingness, you can run this R script but changing the directory path:

```
missing <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/03_indvMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.vcf.imiss", 
                      sep = "\t",
                      header = TRUE)

missing[missing$F_MISS > .2,] # print rows w/ greater than 20% missing
```

Create file to remove the files with missingness >20% called `indv80.txt`

`STL_020` is the only sample with greater than 20% missingness

Also remove `LAE_059` based on high relatedness with other samples.

```
vcftools --vcf ./02_remove_nonbiallelic/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.vcf \
--remove ./03_indvMissing/indv80.txt \
--recode --recode-INFO-all --out ./03_indvMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80
```
***

## Remove sites with > 20% missingness (e.g., 2* 132 idnv *.20 = 52.8, so N = 52)

There are 132 samples (now with the two removed) *.2 * 2 (biallelic loci) so 52 sites missing.

```
vcftools --vcf ./03_indvMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.vcf \
--max-missing-count 52 --recode --recode-INFO-all --out ./04_gntypMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80

```

**Brings number of SNPs to 16,725,129**

***

## Filter  by MAF 

Filter by MAF < 0.05

Filter joint VCF first

```
vcftools --vcf ./04_gntypMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.vcf \
--maf .05 --recode --recode-INFO 'MAF5' \
--out ./05_MAF/joint/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.
```

Loop through file names to do by pop

```
for FILE in LAE LAO STLE STLO3 STLO
do
vcftools --vcf ./04_gntypMissing/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.vcf \
--keep filtered_${FILE}samples.txt \
--maf .05 --recode --recode-INFO-all \
--out ./05_MAF/by_pop/${FILE}samples_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.MAF005
done
```

SNP counts from filtering

| Population      | SNP count |
| ----------- | ----------- |
| LAE     | 5,814,635       |
| LAO  | 6,042,457        |
| STLE     | 4,444,191       |
| STLO3  | 4,253,978        |
| STLO  | 4,552,075        |

## Filter with no MAF but remove invariant sites

See filter_invariantSNPs_byPop.sh file

SNP counts from filtering

| Population  | SNP count   |
| ----------- | ----------- |
| LAE         | 7,167,282   |
| LAO         | 7,337,647   |
| STLE        | 4,996,010   |
| STLO3       | 5,362,161   |
| STLO        | 5,080,474   |

***

## Filter for inversion
vcftools --vcf ../../09_filterVCFs/05_MAF/by_pop/STLOsamples_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.MAF005.recode.vcf --chr NC_060182.1 --from-bp 6403689 --to-bp 35403689 --recode --keep-INFO-all --out STLO_chr10_inversion

vcftools --vcf ../../09_filterVCFs/05_MAF/by_pop/LAOsamples_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.MAF005.recode.vcf --chr NC_060182.1 --from-bp 6403689 --to-bp 35403689 --recode --keep-INFO-all --out LAO_chr10_inversion

