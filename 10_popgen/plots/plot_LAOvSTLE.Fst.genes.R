#=============================================================================================================#
# Script created by Mark Christie, all rights reserved, contact at markchristie1500@gmail.com, ammended by Morgan
# Sparks 2022 (msparks1309@gmail.com)
# Script created in version R 4.0.2 on 12/10/2020
# This script: visualize results of fst analyses
#============================================================================================================#

library(ggplot2); library(patchwork); library(tidyverse); library(ggdist)



setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon")

dat1  <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAOvsSTLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.weir.fst", header = TRUE, sep ="\t") 
head(dat1)
dim(dat1)


#=============================================================================================================#
# Analyses:
# 1. Find Fst outliers with Z scores and no windows



output.full <- dat1

# begin processing for Zfst; remove NAs
length(output.full[, 1])
fst <- output.full[, 3]
length(fst)
fst <- output.full[-(which(is.na(fst) == TRUE)), ]
length(fst[, 1])

# rename scaffolds as a single scaffold name
#scaffolds <- grep("NW_", fst[, 2])
#fst[scaffolds, 2] <- "scaffold"

# plot outliers per chromosome SNP by SNP===============================#
chroms <- unique(fst[, 1])

OUT <- NULL
for(c in 1:length(chroms)){
  
  fst2 <- fst[which(fst[, 1] == chroms[c]), ]
  
  # calculate z score
  fsts     <- fst2$WEIR_AND_COCKERHAM_FST
  mean.fst <- mean(fsts, na.rm = TRUE)
  sd.fst   <- sd(fsts, na.rm = TRUE)
  z.fst <- (fsts-mean.fst)/sd.fst
  
  #plot(fst2$POS, z.fst, main = c)
  output <- cbind(CHROM = fst2$CHROM, POS = fst2$POS, WC_FST = fst2$WEIR_AND_COCKERHAM_FST, z.fst)
  OUT    <- rbind(OUT, output)
  
}

all_snps <- data.frame(OUT)
head(all_snps)

#write.table(all_snps, "./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAOvsSTLO.Zfst.snps.txt", sep = "\t", quote = FALSE)
#write.table(all_snps, "./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAOvsSTLE.Zfst.snps.txt", sep = "\t", quote = FALSE)


# plot outliers per chromosome in sliding windows ===============================#
window.length <- 100000 #100kb
window.step   <- 50000  #50kb

OUT2 <- NULL
for(c in 1:length(chroms)){
  
  fst2 <- fst[which(fst[, 1] == chroms[c]), ]
  
  # calculate z score
  
  window.start <- fst2$POS[1]
  window.end   <- fst2$POS[length(fst2[, 1])]
  
  wstarts <- seq(from = window.start, to = window.end, by = window.step)
  wends   <- wstarts + window.length
  
  OUT3 <- NULL
  for(w in 1:length(wstarts)){
    window   <- wstarts[w]:wends[w]
    fst3     <- fst2[which(fst2$POS >= wstarts[w] & fst2$POS < wends[w]), ]
    fsts     <- fst3$WEIR_AND_COCKERHAM_FST
    mean.fst <- mean(fsts, na.rm = TRUE)
    out      <- cbind(w, wstarts[w], wends[w], mean.fst)
    OUT3     <- rbind(OUT3, out)
  }
  
  # remove nas
  # fsts     <- OUT3[, 4]
  # m1       <- which(is.na(fsts) == TRUE)
  # fsts     <- OUT3[-(which(is.na(fsts) == TRUE)), ]
  
  # remove nas - can occur if snps don't occur in a window
  fsts      <- OUT3
  fsts9     <- fsts[, 4]
  mean.fst  <- mean(fsts9)
  sd.fst    <- sd(fsts9)
  
  z.fst <- (fsts9-mean.fst)/sd.fst
  
  #range(z.fst)
  #plot(fsts[, 2], z.fst, main = c)
  #real.starts <- wstarts[-(m1)]
  real.starts <- wstarts
  ends    <- real.starts + window.length
  output  <- cbind(c,  real.starts, ends, z.fst)
  OUT2    <- rbind(OUT2, output)
  
}


#==========================================================================================================#
# processing and output
head(OUT)
head(OUT2)
output1 <- all_snps[which(as.numeric(all_snps[, 4]) > 5), ]

output2 <- data.frame(OUT2[which(as.numeric(OUT2[, 4]) > 5), ])

zfst <- data.frame(OUT2)


out.windows <- output2

write.table(out.windows,"./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAO_STLE.zfst.windows.txt", sep = "\t")

write.table(output2, "./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAOvSTLEoutlier_windows_all_individuals.txt", sep = "\t")


# plot key outliers -------------------------------------------------------

gff <-  ape::read.gff("./data/assemblies/OgorEven_v1.0/GCF_021184085.1_OgorEven_v1.0_genomic.gff")

gff$seqid <- recode(gff$seqid, 
                    "NC_060173.1" = 1,
                    "NC_060174.1" = 2,
                    "NC_060175.1" = 3,
                    "NC_060176.1" = 4,
                    "NC_060177.1" = 5,
                    "NC_060178.1" = 6,
                    "NC_060179.1" = 7,
                    "NC_060180.1" = 8,
                    "NC_060181.1" = 9,
                    "NC_060182.1" = 10,
                    "NC_060183.1" = 11,
                    "NC_060184.1" = 12,
                    "NC_060185.1" = 13,
                    "NC_060186.1" = 14,
                    "NC_060187.1" = 15,
                    "NC_060188.1" = 16,
                    "NC_060189.1" = 17,
                    "NC_060190.1" = 18,
                    "NC_060191.1" = 19,
                    "NC_060192.1" = 20,
                    "NC_060193.1" = 21,
                    "NC_060194.1" = 22,
                    "NC_060195.1" = 23,
                    "NC_060196.1" = 24,
                    "NC_060197.1" = 25,
                    "NC_060198.1" = 26 )

### window information

outliers.winds.corr <-out.windows

### get Zfst snps for LAO v STLE comparisons

Zfst.snps <- all_snps
Zfst.snps <- Zfst.snps[Zfst.snps$z.fst >= 0,]

### TajD and Het

STLE.tajD <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/TajD/noMAF/STLE_noMAF.TajimaD10KB.Tajima.D", header = TRUE)
STLE.tajD <- chrom.relable(STLE.tajD)

LAO.tajD <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/TajD/noMAF/LAO_noMAF.TajimaD10KB.Tajima.D", header = TRUE)
LAO.tajD <- chrom.relable(LAO.tajD)

#STLO.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/STLO.10KBwindowed.het.txt", header = TRUE)

#STLE.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/", header = TRUE)

LAO.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/LAO.10KBwindowed.het.txt", header = TRUE)


# Plot genes --------------------------------------------------------------

window.1 <- outliers.winds.corr[46,]
snps.1 <-  Zfst.snps[Zfst.snps$CHROM ==window.1$c,]
genes.1 <- gff[gff$seqid ==window.1$c & gff$type %in% c("gene", "CDS"),]
missense.gene.1 <-genes.1[str_which(genes.1$attributes,pattern = "LOC124009961" ),]
missense.gene.1 <- missense.gene.1[missense.gene.1$type == "gene",]
missense.gene.1.1 <-genes.1[str_which(genes.1$attributes,pattern = "lpar5" ),]
missense.gene.1.1 <- missense.gene.1.1[missense.gene.1.1$type == "gene",]
#STLO.het.1 <- STLO.het[STLO.het$chr ==window.1$c,]
LAO.het.1 <- LAO.het[LAO.het$chr ==window.1$c,]
STLE.tajD.1 <- STLE.tajD[STLE.tajD$CHROM == window.1$c,]
LAO.tajD.1 <- LAO.tajD[LAO.tajD$CHROM == window.1$c,]
snps.missense.1 <- snps.1[snps.1$POS %in% c(14729846, 14770211),]
pos.adj.1 <- 400
w.adj <- .1

#per2 mean fst and window fst
mean(snps.1[snps.1$POS >= 14722522 & snps.1$POS <= 14751015, "WC_FST" ])
mean(snps.1$WC_FST)

s.1 <- ggplot() +
  geom_rect(data= missense.gene.1, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_rect(data= missense.gene.1.1, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_point(data = snps.1, aes(x = POS/1e6, y = WC_FST), size = 1, color = "darkgrey") +
  geom_point(data = snps.1[snps.1$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 1) +
  geom_rect(data = snps.missense.1,  aes(xmin = (POS - pos.adj.1)/1e6, xmax= (POS + pos.adj.1)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.1, aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black" )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.1$real.starts/1e6 - w.adj, window.1$ends/1e6 + w.adj)) +
  labs(x = "Chromosome 22 position (Mbp)", y =  bquote(paste(italic('F'['ST']*' '))), subtitle = "period circadian protein homolog 2-like and lysophosphatidic acid receptor 5") +
  lims(y = c(0,1.05)) +
  theme_classic() +
  theme(legend.position = "none")

t.1 <- ggplot() +
  geom_rect(data= missense.gene.1, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLE.tajD.1, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.1, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.1$real.starts/1e6 -w.adj, window.1$ends/1e6 + w.adj)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  bquote(paste("Tajima's ", italic("D")))) +
  theme_classic()

h.1 <- ggplot() +
  geom_rect(data= missense.gene.1, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  #geom_line(data = STLO.het.1, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.1, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.1$real.starts/1e6 -w.adj, window.1$ends/1e6 + w.adj), y = c(0,.5)) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.1$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.9,1),
        legend.background = element_rect(fill = "transparent"))

#### gene 2


window.2 <- outliers.winds.corr[49,]
snps.2 <-  Zfst.snps[Zfst.snps$CHROM ==window.2$c,]
genes.2 <- gff[gff$seqid ==window.2$c & gff$type %in% c("gene", "CDS"),]
missense.gene.2 <- genes.2[str_which(genes.2$attributes,pattern = "LOC124013183" ),]
missense.gene.2 <- missense.gene.2[missense.gene.2$type == "gene",]
#STLO.het.2 <- STLO.het[STLO.het$chr ==window.2$c,]
LAO.het.2 <- LAO.het[LAO.het$chr ==window.2$c,]
STLE.tajD.2 <- STLE.tajD[STLE.tajD$CHROM == window.2$c,]
LAO.tajD.2 <- LAO.tajD[LAO.tajD$CHROM == window.2$c,]
snps.missense.2 <- snps.2[snps.2$POS %in% c(2606585, 2606594, 2607167),]
pos.adj.2 <- 25


s.2 <- ggplot() +
  geom_rect(data= missense.gene.2, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  # geom_rect(data= genes.2[genes.2$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.2$z.fst) +0.5, group = start), alpha = 0.25) +
  geom_point(data = snps.2, aes(x = POS/1e6, y =  WC_FST), size = 1, color = "darkgrey") +
  geom_point(data = snps.2[snps.2$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 1) +
  geom_rect(data = snps.missense.2,   aes(xmin = (POS - pos.adj.2)/1e6, xmax= (POS + pos.adj.2)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.2,aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black", linewidth = .25 )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.2$real.starts/1e6, window.2$ends/1e6)) +
  labs( x = paste("Chromosome", window.2$c, "position (Mbp)"), y =   bquote(paste(italic('F'['ST']*' '))), subtitle = "CD209 antigen-like protein A") +
  theme_classic() +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "transparent"))

# s.2 <- ggplot() +
#   geom_rect(data= genes.2[genes.2$type == "gene",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.2$z.fst) +0.5, group = start), alpha = 0.25) +
#   # geom_rect(data= genes.2[genes.2$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.2$z.fst) +0.5, group = start), alpha = 0.25) +
#   geom_point(data = snps.2, aes(x = POS/1e6, y =  WC_FST), size = 0.25, shape = 1) +
#   geom_point(data = snps.2[snps.2$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 0.25) +
#   geom_rect(data = snps.missense.2,  aes(xmin = (POS - pos.adj.2)/1e6, xmax= (POS + pos.adj.2)/1e6, ymin = max(snps.2$z.fst) , ymax =max(snps.2$z.fst) +0.5), fill = "red") +
#   geom_rect(data = genes.2, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(snps.2$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
#   scale_fill_manual(name = "", values = c("black", "white")) +
#   coord_cartesian( x = c(window.2$starts/1e6+ .075, window.2$ends/1e6 + .03)) +
#   labs(x = NULL, y =  "ZFst", subtitle = "CD209 antigen-like protein A") +
#   theme_classic() +
#   theme(legend.position = "none",
#         legend.background = element_rect(fill = "transparent"))

t.2 <- ggplot() +
  geom_rect(data= missense.gene.2, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLE.tajD.2, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.2, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.2$starts/1e6, window.2$ends/1e6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  "Tajima's D") +
  theme_classic()

h.2 <- ggplot() +
  geom_rect(data= missense.gene.2, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.het.2, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.2, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.2$starts/1e6, window.2$ends/1e6), y = c(0,.5)) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.2$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.1,.9),
        legend.background = element_rect(fill = "transparent"))

### gene 3

window.3 <- outliers.winds.corr[3,]
snps.3 <-  Zfst.snps[Zfst.snps$CHROM ==window.3$c,]
genes.3 <-gff[gff$seqid ==window.3$c & gff$type %in% c("gene", "CDS"),]
missense.gene.3 <- genes.3[str_which(genes.3$attributes,pattern = "gnrhr4" ),]
missense.gene.3 <- missense.gene.3[missense.gene.3$type == "gene",]
#STLO.het.3 <- STLO.het[STLO.het$chr ==window.3$c,]
LAO.het.3 <- LAO.het[LAO.het$chr ==window.3$c,]
STLE.tajD.3 <- STLE.tajD[STLE.tajD$CHROM == window.3$c,]
LAO.tajD.3 <- LAO.tajD[LAO.tajD$CHROM == window.3$c,]
snps.missense.3 <- snps.3[snps.3$POS ==32226620,]
pos.adj.3 <-100

s.3 <- ggplot() +
  geom_rect(data= missense.gene.3, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.35) +
  geom_point(data = snps.3, aes(x = POS/1e6, y =  WC_FST), color = "darkgrey") +
  geom_point(data = snps.3[snps.3$z.fst >= 4,], aes(x = POS/1e6, y = WC_FST)) +
  geom_rect(data = snps.missense.3,  aes(xmin = (POS - pos.adj.1)/1e6, xmax= (POS + pos.adj.1)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.3, aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black", linewidth = .35)+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.3$real.starts/1e6, window.3$ends/1e6)) +
  labs( x = paste("Chromosome", window.3$c, "position (Mbp)"), y =   bquote(paste(italic('F'['ST']*' '))), subtitle = "gonadotropin releasing hormone receptor 4") +
  lims(y = c(0,1.05)) +
  theme_classic() +
  theme(legend.position = "none")

### plot as Zfst
# s.3 <- ggplot() +
#   geom_rect(data= genes.3[genes.3$type == "gene",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.3$z.fst) +0.5, group = start), alpha = 0.35) +
#   # geom_rect(data= genes.3[genes.3$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.3$z.fst) +0.5, group = start), alpha = 0.35) +
#   geom_point(data = snps.3, aes(x = POS/1e6, y =  z.fst), size = 0.35, shape = 1) +
#   geom_point(data = snps.3[snps.3$z.fst >= 4,], aes(x = POS/1e6, y =  z.fst), size = 0.35) +
#   geom_rect(data = snps.missense.3,  aes(xmin = (POS - pos.adj.3)/1e6, xmax= (POS + pos.adj.3)/1e6, ymin = max(snps.3$z.fst) , ymax =max(snps.3$z.fst) +0.5), fill = "red") +
#   geom_rect(data = genes.3, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(snps.3$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
#   scale_fill_manual(name = "", values = c("black", "white")) +
#   coord_cartesian( x = c(window.3$starts/1e6, window.3$ends/1e6)) +
#   labs(x = NULL, y =   bquote(paste('Z',italic('F'['ST']*' '))), subtitle = "gonadotropin releasing hormone receptor 4", title = "c)") +
#   theme_classic() +
#   theme(legend.position = "none")

t.3 <- ggplot() +
  geom_rect(data= missense.gene.3, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.35) +
  geom_line(data = STLE.tajD.3, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.3, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.3$starts/1e6, window.3$ends/1e6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  bquote(paste("Tajima's ", italic("D")))) +
  theme_classic()

h.3 <- ggplot() +
  geom_rect(data= missense.gene.3, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.35) +
  #geom_line(data = STL.het.3, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.3, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.3$starts/1e6, window.3$ends/1e6), y = c(0,.5)) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.3$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.15,1),
        legend.background = element_rect(fill = "transparent"))

ow.3 <- s.3/t.3/h.3

### 

plot.SI <- (s.3/s.1/s.2) + plot_annotation(title = 'BC Odd vs. GL Even')

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/BCOddvGLEven.missensegenes.pdf", plot.SI, height  = 8, width = 8.5, units = "in", dpi = 300)


# Functions ---------------------------------------------------------------

chrom.relable <- function(data){
  
  data$CHROM <-  recode(data$CHROM , 
                        "NC_060173.1" = 1,
                        "NC_060174.1" = 2,
                        "NC_060175.1" = 3,
                        "NC_060176.1" = 4,
                        "NC_060177.1" = 5,
                        "NC_060178.1" = 6,
                        "NC_060179.1" = 7,
                        "NC_060180.1" = 8,
                        "NC_060181.1" = 9,
                        "NC_060182.1" = 10,
                        "NC_060183.1" = 11,
                        "NC_060184.1" = 12,
                        "NC_060185.1" = 13,
                        "NC_060186.1" = 14,
                        "NC_060187.1" = 15,
                        "NC_060188.1" = 16,
                        "NC_060189.1" = 17,
                        "NC_060190.1" = 18,
                        "NC_060191.1" = 19,
                        "NC_060192.1" = 20,
                        "NC_060193.1" = 21,
                        "NC_060194.1" = 22,
                        "NC_060195.1" = 23,
                        "NC_060196.1" = 24,
                        "NC_060197.1" = 25,
                        "NC_060198.1" = 26 )
  return(data)
}

