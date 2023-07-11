# Filtering for run1 samples

**Original unfiltered file has 24,071,073 SNPs**

Query SNP number using `grep -v "^#" FILE.vcf | wc -l &`

## Remove gatk hardfiltered

Take file named GL_PinkSalmon_hardfiltered_snps.vcf.gz and remove gatk flagged hard filtered sites. Out file is `GL_PinkSalmon_hardfiltered_snps.filt.vcf.gz`

```
vcftools --gzvcf ./gatk_HardFilter/run1_allSamps/GL_PinkSalmon_hardfiltered_snps.vcf.gz --remove-filtered-all --recode --recode-INFO-all --out ./soft_filtering/GL_PinkSalmon_hardfiltered_snps.filt.vcf.gz
```

**Number of SNPs = 20,194,631**

## Remove non-biallelic sites

Take file names GL_PinkSalmon_hardfiltered_snps.filt.vcf.gz and remove all non biallelic sites. Out file is 

```
vcftools --gzvcf ./soft_filtering/run1_allSamps/GL_PinkSalmon_hardfiltered_snps.filt.vcf.gz --min-alleles 2 --max-alleles 2  --recode --recode-INFO-all --out ./soft_filtering/run1_allSamps/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.vcf.gz
```

**Number of SNPs = 18,176,215**

## Determine missing idividuals

Make list of genotype missingness by individual


```
vcftools --vcf $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.vcf --missing-indv --out $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.vcf
```


**STL_020 is the only sample with >20% missingness**

Write out file with individual names with 80% genotype presence (determined in R), name it `indiv80.txt`


Remove indviduals with 20% genotype missingness

```
vcftools --vcf $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.vcf \
--keep $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/indiv80.txt \
--recode --recode-INFO-all --out $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.vcf
```

**Number of SNPs = 18176215**


## Remove sites with > 20% missingness (e.g., 2*133*.20 = 53.2, so N = 52)
```
vcftools --vcf $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.vcf \
--missing-site --out $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/gntyp_missingness/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.vcf
```

```
vcftools --vcf $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.vcf \
--max-missing-count 52 --recode --recode-INFO-all --out $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/gntyp_missingness/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.vcf
```
**Number of SNPs = 16,680,075**

## Remove sites with MAF <0.05

```
vcftools --vcf $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/gntyp_missingness/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.vcf \
--maf .05 --recode --recode-INFO 'MAF5' --out $FILTERING/soft_filtering/run1_allSamps/filt_biallelic/indiv_missingness/gntyp_missingness/MAF_filter/GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.
```

**Number of SNPs = 5,331,817**

## Compute post filtering metrics

Compute metrics on individual missingness, genotype missingness, hardy-weinberg eq, allele frequency, and site mean depth
```
vcftools --vcf $FILE --missing-indv --out $FILE
vcftools --vcf $FILE --missing-site --out $FILE
vcftools --vcf $FILE --hardy --out $FILE
vcftools --vcf $FILE --freq --out $FILE
vcftools --vcf $FILE --site-mean-depth --out $FILE
```

## Split into populations specific vcfs to filter MAF by population
```
vcftools --vcf ../GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.vcf --maf .05 --keep filtered_STLO3samples.txt --recode --recode-INFO-all --out STLO3samples_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf
```
**Number of SNPs = 4,258,789**

```
vcftools --vcf ../GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.vcf --maf .05 --keep filtered_STLEsamples.txt --recode --recode-INFO-all --out STLEsamples_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf
```
**Number of SNPs = 4,450,204**

```
vcftools --vcf ../GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.vcf --maf .05 --keep filtered_LAOsamples.txt --recode --recode-INFO-all --out LAOsamples_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf
```
**Number of SNPs = 6,031,542**

```
vcftools --vcf ../GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.vcf --maf .05 --keep filtered_LAEsamples.txt --recode --recode-INFO-all --out LAEsamples_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf
```
**Number of SNPs = 6,077,481**

```
vcftools --vcf ../GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.vcf --maf .05 --keep filtered_STLOsamples.txt --recode --recode-INFO-all --out STLOsamples_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.vcf
```
**Number of SNPs = 4,555,518**