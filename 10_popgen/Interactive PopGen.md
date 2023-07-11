# Interactive PopGen

### Calculate heterozygosity

Will calculate from vcf uploaded into R with `gaston` package. 

First Navigate to folder with vcf.

```
library(gaston); library(tidyverse)

# read in vcf
GL.vcf <- read.vcf("./GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.MAF005.recode.vcf", )

# create bed matrix
GL.bed <- as.matrix(GL.vcf)

# make into a df and separate Samples into pops and sample #
GL.bed <- data.frame(cbind(samples = rownames(GL.bed), GL.bed))
GL.bed <- GL.bed %>% separate(samples, into = c("Pop", "Sample"), sep = "_")

# function to calculate by pop het across loci

by.locus.het <- function(x){
    sub.df <- as.matrix(GL.bed[GL.bed$Pop == x,3:ncol(GL.bed)]) # subset on pop
    het <- sapply(1:ncol(sub.df), FUN = function(x) length(sub.df[sub.df[,x]==1, x])/length(sub.df[is.na(sub.df[,x]) ==FALSE, x])) # use sapply to calculate per locus het (eg. # of 1s over total number of genotypes that aren't NA)
    return(mean(het))
}

# make an object with het returned
GL.bylocus.het <- sapply(unique(GL.bed$Pop), FUN = by.locus.het)

### calcukate across individuals (alot easier with the ped object from gaston)

GL.byindv.het <- GL.vcf@ped %>%
  separate(samples, into = c("Pop", "Sample") %>% 
             group_by(Pop) %>% 
             summarise(het = mean(hz))

```
