#=============================================================================================================#
# Script created by Mark Christie, all rights reserved, contact at markchristie1500@gmail.com
# Script created in version R 4.1.2 on 2/8/2022
# This script: takes a vcf file and calculates HWE p-values per population
# Usage notes: 
#============================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd(".")
library('gaston') # had to first install rcppp via command line
library('HardyWeinberg') 

list.files()

dat  <- read.vcf("GL_PinkSalmon_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.newchroms.vcf", convert.chr = FALSE) # this file has all the correct indivdiual ids and population assignments
dat
dim(dat)

#write.table(offs, "temp".txt", col.names = TRUE, sep="\t", append = FALSE)
#=============================================================================================================#

# Start isolating populations
dat1 <- as.matrix(dat)
ids  <- rownames(dat1)

pops <- unlist(strsplit(ids, "_"))
pops <- pops[seq(from = 1, to = length(pops), by = 2)]  # formatting here will be very VCF file specific
table(pops)

# create unique populations
popids      <- pops # real pops with yoy separated is included in info
unique.pops <- unique(popids)


OUT <- NULL
for(n in 1:length(unique.pops)){
  upop   <- unique.pops[n]
  tpop   <- dat1[popids == upop, ]
  
  #isolate individual loci
  nloci <- ncol(tpop)
  OUT2 <- NULL
  for(l in 1:nloci){
    l1   <- tpop[, l]
    l2   <- data.frame(table(l1))
    l2   <- l2[, 2]           # l2 must have observations in the order (a_11, a_21, a_22, a_31, ..., a_kk\) 
    if(length(l2) < 2) {next} # if locus is monomorphic (which can happen in individual pops), skip
    if(length(l2)  == 2) {    # format biallelic loci with only 2 genotypes to have 3 genotypes (add 0 to appropriate position) 
      l2     <- data.frame(table(l1))
      l2temp <- cbind(c(0,1,2), c(0,0,0))
      m1     <- match(l2[, 1], l2temp)
      l2temp[m1, 2] <- l2[, 2]
      l2     <- l2temp[, 2]             
    }       
    test <- HWExact(l2, verbose = FALSE)
    out <- cbind(l, test$pval)
    OUT2 <- rbind(OUT2, out)    
  }
  out    <- cbind(upop, OUT2)
  OUT    <- rbind(OUT, out)
}


# Find loci that are out of HWE
hwe   <- OUT
alpha <- 0.05/ncol(dat1) # total number of loci; may need to adjust if lots of monomorphic loci

hwe2  <- hwe[which(as.numeric(hwe[, 3]) < alpha), ] # which values less than alpha
table(hwe2[, 1])
((table(hwe2[, 1]))/ncol(dat1)) * 100 # percent of loci out of HWE

# range is 0.01% to 0.89% (all less than 1%!)


out_per_pop <- data.frame(table(hwe2[, 1]))

out_per_loc <- data.frame(table(hwe2[, 2]))
out_per_loc <-out_per_loc[order(out_per_loc[, 2], decreasing = TRUE), ]

# 1. Average number of loci out of HWE per population
avg <- mean(out_per_pop[, 2])
avg
avg/ncol(dat1)

# STOPPED HERE FOR YPERCH RAD SEQ
# BELOW CODE COULD BE USED TO REMOVE HWE DEVIATES == NOT NEEDED HERE

# 2. Total number of loci out of HWE in 13/26 

#nloc <- which(out_per_loc[, 2] > 18)
##nloc <- length(nloc)
#nloc
#nloc/3312
#18/26

#3. export vcf file with all loci removed that occurred in greater than 3 populations (e.g, 4 or more)
nloc <- which(out_per_loc[, 2] > 3)
loci <- out_per_loc[nloc, 1]

loci <- as.numeric(as.character(loci))
dat2 <- dat1[, -loci]


write.table(dat2, "hwe.vcf", row.names = TRUE, col.names = TRUE, sep="\t", append = FALSE)

loci <- cbind(loci, colnames(dat1)[loci])
write.table(loci, "./HWE_excess/hwe_loci.txt", row.names = TRUE, col.names = TRUE, sep="\t", append = FALSE)


# colnames(dat1)[loci]
#                ]