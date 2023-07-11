library(gdsfmt); library(SNPRelate)

### read in files
vcf <- "/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/05_MAF/joint/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.MAF005.recode.newnames.newchr.vcf"

snpgdsVCF2GDS(vcf, "pink_salmon.joint.filtered.gds", method="biallelic.only")

genofile <- snpgdsOpen("pink_salmon.joint.filtered.gds")

samp.id <- seqGetData(genofile, "sample.id")

pop_code <- scan("/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/PCA/SNPrelate.samples.txt", what=character())

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

### prune snps and run Fst

set.seed(1000)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)

# 123,924 markers are selected in total.

#compute pairwise Fst
pops <- unique(pop_code)

OUT =NULL
for(p1 in pops){
    OUT2 =NULL
    pops2 <- pops[pops != p1]
    for (p2 in pops2){
        flag <- pop_code %in% c(p1, p2)
        samp.sel <- sample.id[flag]
        pop.sel <- pop_code[flag]
        v <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), method="W&C84")
        temp.row <- cbind(p1,p2, v$MeanFst)
        OUT2 <- rbind(OUT2, temp.row)
    }
    OUT <- rbind(OUT, OUT2)
}

pairwise.fst <- data.frame(OUT)
colnames(pairwise.fst) <- c("pop1", "pop2", "fst")
pairwise.fst$fst <- as.numeric(as.character(pairwise.fst$fst))

# get unique values
pairwise.fst.cor <- pairwise.fst[!duplicated(pairwise.fst$fst), ]
pairwise.fst.cor$fst <- round(pairwise.fst.cor$fst, 5)
