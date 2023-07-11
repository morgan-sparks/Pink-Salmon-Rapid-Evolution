# Header ------------------------------------------------------------------
### libraries
library(tidyverse)
###setwd

setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon")

### read gff and change chrom names

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

outliers.winds.corr <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers/LAOvSTLO_outlierwindows_corrected.txt", header = TRUE)


### get Zfst snps for LAO v STLO comparisons

Zfst.snps <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAOvsSTLO.Zfst.snps.txt", sep = "\t", header = TRUE)

Zfst.snps <- data.frame(Zfst.snps)
Zfst.snps <- Zfst.snps[Zfst.snps$z.fst >= 0,]

### Get raw Fst

raw.fst <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAOvsSTLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.weir.fst", sep = "\t", header = TRUE)


### TajD and Het

STLO.tajD <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/TajD/noMAF/STLO_noMAF.TajimaD10KB.Tajima.D", header = TRUE)
STLO.tajD <- chrom.relable(STLO.tajD)

LAO.tajD <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/TajD/noMAF/LAO_noMAF.TajimaD10KB.Tajima.D", header = TRUE)
LAO.tajD <- chrom.relable(LAO.tajD)

STLO.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/STLO.10KBwindowed.het.txt", header = TRUE)

LAO.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/LAO.10KBwindowed.het.txt", header = TRUE)

# Get best outliers -------------------------------------------------------

 ZFst.outs.sum <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers/STLO_ZFstSNP.inwindows.snpEff_genes.txt", sep = "\t", header = TRUE)

 ZFst.outs.effect <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers/snpEff_results/STLO_snpEff_ZfstSNPs_inZfstwindows.txt", sep = "\t", header = TRUE) 
 
Zfst.missense <- ZFst.outs.effect %>% 
  filter(ANN....EFFECT == "missense_variant") %>%
  group_by(EFF....GENE) %>% 
  tally() %>% 
  arrange(desc(n))

Zfst.missense.all <- ZFst.outs.effect %>% 
  filter(ANN....EFFECT == "missense_variant") %>%
  group_by(EFF....GENE) %>% 
  distinct(POS, .keep_all = TRUE)

Zfst.per2 <- ZFst.outs.effect %>% 
  filter(EFF....GENE == "LOC124009961") %>%
  distinct(POS, .keep_all = TRUE)

### get gene posiitions

OUT.genes <- NULL
for (g in 1:nrow(outliers.winds.corr)){

  window <- outliers.winds.corr[g,]
  window.s <- window$starts
  window.e <- window$ends
  window.c <- window$c
  
  inside <- gff[gff$seqid == window.c & gff$start >= window.s & gff$end <= window.e,]
  
  start.overlap <- gff[gff$seqid == window.c & gff$start <= window.s & gff$end >= window.s,]
  
  end.overlap <-   gff[gff$seqid == window.c & gff$start <= window.e & gff$end >= window.e,]

  genes <- rbind(inside, start.overlap, end.overlap)
  
  OUT.genes <- rbind(OUT.genes, genes)
}

Zfst.gene.pos <- data.frame(OUT.genes)

Zfst.gene.pos.coding <-  Zfst.gene.pos %>% 
  filter(type %in% c("gene", "CDS"))


# Plot windows/genes --------------------------------------------------------------

### window.1 
window.1 <- outliers.winds.corr[31,]
snps.1 <-  Zfst.snps[Zfst.snps$CHROM ==window.1$c,]
genes.1 <- gff[gff$seqid ==window.1$c & gff$type %in% c("gene", "CDS"),]
missense.gene.1 <- genes.1[str_which(genes.1$attributes,pattern = "LOC124009961" ),]
missense.gene.1 <- missense.gene.1[missense.gene.1$type == "gene",]
STLO.het.1 <- STLO.het[STLO.het$chr ==window.1$c,]
LAO.het.1 <- LAO.het[LAO.het$chr ==window.1$c,]
STLO.tajD.1 <- STLO.tajD[STLO.tajD$CHROM == window.1$c,]
LAO.tajD.1 <- LAO.tajD[LAO.tajD$CHROM == window.1$c,]
snps.missense.1 <- snps.1[snps.1$POS %in% c(14729846, 14770211),]
pos.adj.1 <- 400
w.adj <- .1

#per2 mean fst and window fst
mean(snps.1[snps.1$POS >= 14722522 & snps.1$POS <= 14751015, "WC_FST" ])
mean(snps.1$WC_FST)

s.1 <- ggplot() +
  geom_rect(data= missense.gene.1, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_point(data = snps.1, aes(x = POS/1e6, y = WC_FST), color = "darkgrey", size = 0.5) +
  geom_point(data = snps.1[snps.1$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 0.5) +
  geom_rect(data = snps.missense.1,  aes(xmin = (POS - pos.adj.1)/1e6, xmax= (POS + pos.adj.1)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.1, aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black", linewidth = 0.25)+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.1$starts, window.1$ends)/1e6) +
  labs(x = NULL, y =  bquote(paste(italic('F'['ST']*' '))), subtitle = "period circadian protein homolog 2-like", title = "d)") +
  lims(y = c(0,1.05)) +
  theme_classic() +
  theme(legend.position = "none")


### plot as Zfst
# ggplot() +
#   geom_rect(data= genes.1[genes.1$type == "gene",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.1$z.fst) +0.5, group = start), alpha = 0.25) +
#   #geom_rect(data= genes.1[genes.1$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.1$z.fst) +0.5, group = start), alpha = 0.25) +
#   geom_point(data = snps.1, aes(x = POS/1e6, y =  z.fst), size = 1, shape = 1) +
#   geom_point(data = snps.1[snps.1$z.fst >= 4,], aes(x = POS/1e6, y =  z.fst), size = 1) +
#   geom_rect(data = snps.missense.1,  aes(xmin = (POS - pos.adj.1)/1e6, xmax= (POS + pos.adj.1)/1e6, ymin = max(snps.1$z.fst) , ymax =max(snps.1$z.fst) +0.5), fill = "red") +
#   geom_rect(data = genes.1, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(snps.1$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
#   scale_fill_manual(name = "", values = c("black", "white")) +
#   coord_cartesian( x = c(window.1$starts/1e6, window.1$ends/1e6)) +
#   labs(x = NULL, y =  bquote(paste('Z',italic('F'['ST']*' '))), subtitle = "period circadian protein homolog 2-like", title = "d)") +
#   theme_classic() +
#   theme(legend.position = "none")

t.1 <- ggplot() +
  geom_rect(data= missense.gene.1, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.tajD.1, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.1, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.1$starts, window.1$ends)/1e6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  bquote(paste("Tajima's ", italic("D")))) +
  theme_classic()

h.1 <- ggplot() +
  geom_rect(data= missense.gene.1, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.het.1, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.1, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.1$starts, window.1$ends)/1e6) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.1$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.9,1),
        legend.background = element_rect(fill = "transparent"))

ow.1 <- s.1/t.1/h.1

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/Fig.4.d.pdf", ow.1, height = 5, width = 3.75, units = "in", dpi  = 600)


### window.2

window.2 <- outliers.winds.corr[3,]
snps.2 <-  Zfst.snps[Zfst.snps$CHROM ==window.2$c,]
genes.2 <- gff[gff$seqid ==window.2$c & gff$type %in% c("gene", "CDS"),]
missense.gene.2 <- genes.2[str_which(genes.2$attributes,pattern = "gnrhr4" ),]
missense.gene.2 <- missense.gene.2[missense.gene.2$type == "gene",]
STLO.het.2 <- STLO.het[STLO.het$chr ==window.2$c,]
LAO.het.2 <- LAO.het[LAO.het$chr ==window.2$c,]
STLO.tajD.2 <- STLO.tajD[STLO.tajD$CHROM == window.2$c,]
LAO.tajD.2 <- LAO.tajD[LAO.tajD$CHROM == window.2$c,]
snps.missense.2 <- snps.2[snps.2$POS ==32226620,]
pos.adj.2 <-50

s.2 <- ggplot() +
  geom_rect(data= missense.gene.2, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_point(data = snps.2, aes(x = POS/1e6, y =  WC_FST), color = "darkgrey", size = 0.5) +
  geom_point(data = snps.2[snps.2$z.fst >= 4,], aes(x = POS/1e6, y = WC_FST), size = 0.5) +
  geom_rect(data = snps.missense.2,  aes(xmin = (POS - pos.adj.2)/1e6, xmax= (POS + pos.adj.2)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.2, aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black", linewidth = .25)+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.2$starts, window.2$ends)/1e6) +
  labs( x = paste("Chromosome", window.2$c, "position (Mbp)"), y =   bquote(paste(italic('F'['ST']*' '))), subtitle = "gonadotropin releasing hormone receptor 4") +
  lims(y = c(0,1.05)) +
  theme_classic() +
  theme(legend.position = "none")

### plot as Zfst
# s.2 <- ggplot() +
#   geom_rect(data= genes.2[genes.2$type == "gene",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.2$z.fst) +0.5, group = start), alpha = 0.25) +
#   # geom_rect(data= genes.2[genes.2$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.2$z.fst) +0.5, group = start), alpha = 0.25) +
#   geom_point(data = snps.2, aes(x = POS/1e6, y =  z.fst), size = 0.25, shape = 1) +
#   geom_point(data = snps.2[snps.2$z.fst >= 4,], aes(x = POS/1e6, y =  z.fst), size = 0.25) +
#   geom_rect(data = snps.missense.2,  aes(xmin = (POS - pos.adj.2)/1e6, xmax= (POS + pos.adj.2)/1e6, ymin = max(snps.2$z.fst) , ymax =max(snps.2$z.fst) +0.5), fill = "red") +
#   geom_rect(data = genes.2, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(snps.2$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
#   scale_fill_manual(name = "", values = c("black", "white")) +
#   coord_cartesian( x = c(window.2$starts/1e6, window.2$ends/1e6)) +
#   labs(x = NULL, y =   bquote(paste('Z',italic('F'['ST']*' '))), subtitle = "gonadotropin releasing hormone receptor 4", title = "c)") +
#   theme_classic() +
#   theme(legend.position = "none")

t.2 <- ggplot() +
  geom_rect(data= missense.gene.2, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.tajD.2, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.2, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.2$starts-5e4, window.2$ends+5e4)/1e6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  bquote(paste("Tajima's ", italic("D")))) +
  theme_classic()

h.2 <- ggplot() +
  geom_rect(data= missense.gene.2, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.het.2, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.2, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.2$starts-5e4, window.2$ends+5e4)/1e6) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.2$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.15,1),
        legend.background = element_rect(fill = "transparent"))

ow.2 <- s.2/t.2/h.2

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/Fig.4.c.pdf", ow.2, height = 5, width = 3.75, units = "in", dpi  = 600)


### window.3

window.3 <- outliers.winds.corr[32,]
snps.3 <-  Zfst.snps[Zfst.snps$CHROM ==window.3$c,]
genes.3 <- gff[gff$seqid ==window.3$c & gff$type %in% c("gene", "CDS"),]
missense.gene.3 <- genes.3[str_which(genes.3$attributes,pattern = "LOC124013183" ),]
missense.gene.3 <- missense.gene.3[missense.gene.3$type == "gene",]
STLO.het.3 <- STLO.het[STLO.het$chr ==window.3$c,]
LAO.het.3 <- LAO.het[LAO.het$chr ==window.3$c,]
STLO.tajD.3 <- STLO.tajD[STLO.tajD$CHROM == window.3$c,]
LAO.tajD.3 <- LAO.tajD[LAO.tajD$CHROM == window.3$c,]
snps.missense.3 <- snps.3[snps.3$POS %in% c(2606585, 2606594, 2607167),]
pos.adj.3 <- 100


s.3 <- ggplot() +
  geom_rect(data= missense.gene.3, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  # geom_rect(data= genes.3[genes.3$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.3$z.fst) +0.5, group = start), alpha = 0.35) +
  geom_point(data = snps.3, aes(x = POS/1e6, y =  WC_FST), color = "darkgrey", size =0.5) +
  geom_point(data = snps.3[snps.3$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size =0.5) +
  geom_rect(data = snps.missense.3,   aes(xmin = (POS - pos.adj.3)/1e6, xmax= (POS + pos.adj.3)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.3,aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black", linewidth = .25 )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.3$starts-5e4, window.3$ends+5e4)/1e6) +
  lims(y = c(0,1.05))+
  labs( x = paste("Chromosome", window.3$c, "position (Mbp)"), y =   bquote(paste(italic('F'['ST']*' '))), subtitle = "CD209 antigen-like protein A") +
  theme_classic() +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "transparent"))

# s.3 <- ggplot() +
#   geom_rect(data= genes.3[genes.3$type == "gene",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.3$z.fst) +0.5, group = start), alpha = 0.35) +
#   # geom_rect(data= genes.3[genes.3$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.3$z.fst) +0.5, group = start), alpha = 0.35) +
#   geom_point(data = snps.3, aes(x = POS/1e6, y =  WC_FST), size = 0.25, shape = 1) +
#   geom_point(data = snps.3[snps.3$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 0.25) +
#   geom_rect(data = snps.missense.3,  aes(xmin = (POS - pos.adj.3)/1e6, xmax= (POS + pos.adj.3)/1e6, ymin = max(snps.3$z.fst) , ymax =max(snps.3$z.fst) +0.5), fill = "red") +
#   geom_rect(data = genes.3, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(snps.3$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
#   scale_fill_manual(name = "", values = c("black", "white")) +
#   coord_cartesian( x = c(window.3$starts/1e6+ .075, window.3$ends/1e6 + .03)) +
#   labs(x = NULL, y =  "ZFst", subtitle = "CD209 antigen-like protein A") +
#   theme_classic() +
#   theme(legend.position = "none",
#         legend.background = element_rect(fill = "transparent"))

t.3 <- ggplot() +
  geom_rect(data= missense.gene.3, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.tajD.3, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.3, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.3$starts-5e4, window.3$ends+5e4)/1e6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  "Tajima's D") +
  theme_classic()

h.3 <- ggplot() +
  geom_rect(data= missense.gene.3, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.het.3, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.3, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.3$starts-5e4, window.3$ends+5e4)/1e6) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.3$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.1,.9),
        legend.background = element_rect(fill = "transparent"))

ow.3 <- s.3/t.3/h.3


ow.2 | ow.1

### window 4

window.4 <- outliers.winds.corr[11,]
snps.4 <-  Zfst.snps[Zfst.snps$CHROM ==window.4$c,]
genes.4 <- gff[gff$seqid ==window.4$c & gff$type %in% c("gene", "CDS"),]
missense.gene.4 <- genes.4[str_which(genes.4$attributes,pattern = "LOC124043346" ),]
missense.gene.4 <- missense.gene.4[missense.gene.4$type == "gene",]
STLO.het.4 <- STLO.het[STLO.het$chr ==window.4$c,]
LAO.het.4 <- LAO.het[LAO.het$chr ==window.4$c,]
STLO.tajD.4 <- STLO.tajD[STLO.tajD$CHROM == window.4$c,]
LAO.tajD.4 <- LAO.tajD[LAO.tajD$CHROM == window.4$c,]
snps.missense.4 <- snps.4[snps.4$POS %in% c(28065291),]
pos.adj.4 <- 300


s.4 <- ggplot() +
  geom_rect(data= missense.gene.4, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  # geom_rect(data= genes.4[genes.4$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.4$z.fst) +0.5, group = start), alpha = 0.45) +
  geom_point(data = snps.4, aes(x = POS/1e6, y =  WC_FST), size = 1, color = "darkgrey") +
  geom_point(data = snps.4[snps.4$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 1) +
  geom_rect(data = snps.missense.4,   aes(xmin = (POS - pos.adj.4)/1e6, xmax= (POS + pos.adj.4)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.4,aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black", linewidth = .25 )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.4$starts/1e6, window.4$ends/1e6)) +
  labs(x = paste("Chromosome", window.4$c, "position (Mbp)"), y =   bquote(paste(italic('F'['ST']*' '))), subtitle = "collagen alpha-2(V) chain-like") +
  theme_classic() +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "transparent"))

# s.4 <- ggplot() +
#   geom_rect(data= genes.4[genes.4$type == "gene",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.4$z.fst) +0.5, group = start), alpha = 0.45) +
#   # geom_rect(data= genes.4[genes.4$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.4$z.fst) +0.5, group = start), alpha = 0.45) +
#   geom_point(data = snps.4, aes(x = POS/1e6, y =  WC_FST), size = 0.25, shape = 1) +
#   geom_point(data = snps.4[snps.4$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 0.25) +
#   geom_rect(data = snps.missense.4,  aes(xmin = (POS - pos.adj.4)/1e6, xmax= (POS + pos.adj.4)/1e6, ymin = max(snps.4$z.fst) , ymax =max(snps.4$z.fst) +0.5), fill = "red") +
#   geom_rect(data = genes.4, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(snps.4$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
#   scale_fill_manual(name = "", values = c("black", "white")) +
#   coord_cartesian( x = c(window.4$starts/1e6+ .075, window.4$ends/1e6 + .03)) +
#   labs(x = NULL, y =  "ZFst", subtitle = "CD209 antigen-like protein A") +
#   theme_classic() +
#   theme(legend.position = "none",
#         legend.background = element_rect(fill = "transparent"))

t.4 <- ggplot() +
  geom_rect(data= missense.gene.4, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.tajD.4, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.4, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.4$starts/1e6, window.4$ends/1e6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  "Tajima's D") +
  theme_classic()

h.4 <- ggplot() +
  geom_rect(data= missense.gene.4, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.het.4, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.4, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.4$starts/1e6, window.4$ends/1e6), y = c(0,.5)) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.4$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.1,.9),
        legend.background = element_rect(fill = "transparent"))

ow.4 <- s.4/t.4/h.4

### window 5

window.5 <- outliers.winds.corr[22,]
snps.5 <-  Zfst.snps[Zfst.snps$CHROM ==window.5$c,]
genes.5 <- gff[gff$seqid ==window.5$c & gff$type %in% c("gene", "CDS"),]
missense.gene.5 <- genes.5[str_which(genes.5$attributes,pattern = "LOC124000748" ),]
missense.gene.5 <- missense.gene.5[missense.gene.5$type == "gene",]
STLO.het.5 <- STLO.het[STLO.het$chr ==window.5$c,]
LAO.het.5 <- LAO.het[LAO.het$chr ==window.5$c,]
STLO.tajD.5 <- STLO.tajD[STLO.tajD$CHROM == window.5$c,]
LAO.tajD.5 <- LAO.tajD[LAO.tajD$CHROM == window.5$c,]
snps.missense.5 <- snps.5[snps.5$POS %in% c(80387446),]
pos.adj.5 <- 150


s.5 <- ggplot() +
  geom_rect(data= missense.gene.5, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  # geom_rect(data= genes.5[genes.5$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.5$z.fst) +0.5, group = start), alpha = 0.55) +
  geom_point(data = snps.5, aes(x = POS/1e6, y =  WC_FST), size = 1, color = "darkgrey") +
  geom_point(data = snps.5[snps.5$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size =1) +
  geom_rect(data = snps.missense.5,   aes(xmin = (POS - pos.adj.5)/1e6, xmax= (POS + pos.adj.5)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.5,aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black", linewidth = .25 )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.5$starts/1e6, window.5$ends/1e6)) +
  labs(x = paste("Chromosome", window.5$c, "position (Mbp)"), y =   bquote(paste(italic('F'['ST']*' '))), subtitle = "B-cell CLL/lymphoma 7 protein family member A-like") +
  theme_classic() +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "transparent"))

# s.5 <- ggplot() +
#   geom_rect(data= genes.5[genes.5$type == "gene",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.5$z.fst) +0.5, group = start), alpha = 0.55) +
#   # geom_rect(data= genes.5[genes.5$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.5$z.fst) +0.5, group = start), alpha = 0.55) +
#   geom_point(data = snps.5, aes(x = POS/1e6, y =  WC_FST), size = 0.25, shape = 1) +
#   geom_point(data = snps.5[snps.5$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 0.25) +
#   geom_rect(data = snps.missense.5,  aes(xmin = (POS - pos.adj.5)/1e6, xmax= (POS + pos.adj.5)/1e6, ymin = max(snps.5$z.fst) , ymax =max(snps.5$z.fst) +0.5), fill = "red") +
#   geom_rect(data = genes.5, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(snps.5$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
#   scale_fill_manual(name = "", values = c("black", "white")) +
#   coord_cartesian( x = c(window.5$starts/1e6+ .075, window.5$ends/1e6 + .03)) +
#   labs(x = NULL, y =  "ZFst", subtitle = "CD209 antigen-like protein A") +
#   theme_classic() +
#   theme(legend.position = "none",
#         legend.background = element_rect(fill = "transparent"))

t.5 <- ggplot() +
  geom_rect(data= missense.gene.5, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.tajD.5, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.5, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.5$starts/1e6, window.5$ends/1e6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  "Tajima's D") +
  theme_classic()

h.5 <- ggplot() +
  geom_rect(data= missense.gene.5, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.het.5, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.5, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.5$starts/1e6, window.5$ends/1e6), y = c(0,.5)) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.5$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.1,.9),
        legend.background = element_rect(fill = "transparent"))

ow.5 <- s.5/t.5/h.5


### window 6

window.6 <- outliers.winds.corr[26,]
snps.6 <-  Zfst.snps[Zfst.snps$CHROM ==window.6$c,]
genes.6 <- gff[gff$seqid ==window.6$c & gff$type %in% c("gene", "CDS"),]
missense.gene.6 <- genes.6[str_which(genes.6$attributes,pattern = "LOC124004644" ),]
missense.gene.6 <- missense.gene.6[missense.gene.6$type == "gene",]
STLO.het.6 <- STLO.het[STLO.het$chr ==window.6$c,]
LAO.het.6 <- LAO.het[LAO.het$chr ==window.6$c,]
STLO.tajD.6 <- STLO.tajD[STLO.tajD$CHROM == window.6$c,]
LAO.tajD.6 <- LAO.tajD[LAO.tajD$CHROM == window.6$c,]
snps.missense.6 <- snps.6[snps.6$POS %in% c(4436398),]
pos.adj.6 <- 300


s.6 <- ggplot() +
  geom_rect(data= missense.gene.6, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  # geom_rect(data= genes.6[genes.6$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.6$z.fst) +0.6, group = start), alpha = 0.65) +
  geom_point(data = snps.6, aes(x = POS/1e6, y =  WC_FST), size = 1, color = "darkgrey") +
  geom_point(data = snps.6[snps.6$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 1) +
  geom_rect(data = snps.missense.6,   aes(xmin = (POS - pos.adj.6)/1e6, xmax= (POS + pos.adj.6)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.6,aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black", linewidth = .25 )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.6$starts/1e6, window.6$ends/1e6)) +
  labs(x = paste("Chromosome", window.6$c, "position (Mbp)"), y =   bquote(paste(italic('F'['ST']*' '))), subtitle = "cytochrome P450 11B, mitochondrial-like") +
  theme_classic() +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "transparent"))

# s.6 <- ggplot() +
#   geom_rect(data= genes.6[genes.6$type == "gene",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.6$z.fst) +0.6, group = start), alpha = 0.65) +
#   # geom_rect(data= genes.6[genes.6$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.6$z.fst) +0.6, group = start), alpha = 0.65) +
#   geom_point(data = snps.6, aes(x = POS/1e6, y =  WC_FST), size = 0.25, shape = 1) +
#   geom_point(data = snps.6[snps.6$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 0.25) +
#   geom_rect(data = snps.missense.6,  aes(xmin = (POS - pos.adj.6)/1e6, xmax= (POS + pos.adj.6)/1e6, ymin = max(snps.6$z.fst) , ymax =max(snps.6$z.fst) +0.6), fill = "red") +
#   geom_rect(data = genes.6, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(snps.6$z.fst) +0.6 , ymax= Inf, fill= type), color = "black" )+
#   scale_fill_manual(name = "", values = c("black", "white")) +
#   coord_cartesian( x = c(window.6$starts/1e6+ .075, window.6$ends/1e6 + .03)) +
#   labs(x = NULL, y =  "ZFst", subtitle = "CD209 antigen-like protein A") +
#   theme_classic() +
#   theme(legend.position = "none",
#         legend.background = element_rect(fill = "transparent"))

t.6 <- ggplot() +
  geom_rect(data= missense.gene.6, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.tajD.6, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.6, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.6$starts/1e6, window.6$ends/1e6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  "Tajima's D") +
  theme_classic()

h.6 <- ggplot() +
  geom_rect(data= missense.gene.6, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.het.6, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.6, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.6$starts/1e6, window.6$ends/1e6), y = c(0,.6)) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.6$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.1,.9),
        legend.background = element_rect(fill = "transparent"))

ow.6 <- s.6/t.6/h.6

### window 7
window.7 <- outliers.winds.corr[31,]
snps.7 <-  Zfst.snps[Zfst.snps$CHROM ==window.7$c,]
genes.7 <- gff[gff$seqid ==window.7$c & gff$type %in% c("gene", "CDS"),]
missense.gene.7 <- genes.7[str_which(genes.7$attributes,pattern = "lpar5a" ),]
missense.gene.7 <- missense.gene.7[missense.gene.7$type == "gene",]
STLO.het.7 <- STLO.het[STLO.het$chr ==window.7$c,]
LAO.het.7 <- LAO.het[LAO.het$chr ==window.7$c,]
STLO.tajD.7 <- STLO.tajD[STLO.tajD$CHROM == window.7$c,]
LAO.tajD.7 <- LAO.tajD[LAO.tajD$CHROM == window.7$c,]
snps.missense.7 <- snps.7[snps.7$POS %in% c(14770211),]
pos.adj.7 <- 300


s.7 <- ggplot() +
  geom_rect(data= missense.gene.7, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  # geom_rect(data= genes.7[genes.7$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.7$z.fst) +0.7, group = start), alpha = 0.75) +
  geom_point(data = snps.7, aes(x = POS/1e6, y =  WC_FST), size = 1, color = "darkgrey") +
  geom_point(data = snps.7[snps.7$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 1) +
  geom_rect(data = snps.missense.7,   aes(xmin = (POS - pos.adj.7)/1e6, xmax= (POS + pos.adj.7)/1e6, ymin = .9, ymax =1.05), fill = "red") +
  geom_rect(data = genes.7,aes(xmin = start/1e6, xmax = end/1e6, ymin = 1.05, ymax =Inf, fill= type), color = "black", linewidth = .25 )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  coord_cartesian( x = c(window.7$starts/1e6, window.7$ends/1e6)) +
  labs(x = paste("Chromosome", window.7$c, "position (Mbp)"), y =   bquote(paste(italic('F'['ST']*' '))), subtitle = "cysophosphatidic acid receptor 5a") +
  theme_classic() +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "transparent"))

# s.7 <- ggplot() +
#   geom_rect(data= genes.7[genes.7$type == "gene",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.7$z.fst) +0.7, group = start), alpha = 0.75) +
#   # geom_rect(data= genes.7[genes.7$type == "CDS",], aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = max(snps.7$z.fst) +0.7, group = start), alpha = 0.75) +
#   geom_point(data = snps.7, aes(x = POS/1e6, y =  WC_FST), size = 0.25, shape = 1) +
#   geom_point(data = snps.7[snps.7$z.fst >= 4,], aes(x = POS/1e6, y =  WC_FST), size = 0.25) +
#   geom_rect(data = snps.missense.7,  aes(xmin = (POS - pos.adj.7)/1e6, xmax= (POS + pos.adj.7)/1e6, ymin = max(snps.7$z.fst) , ymax =max(snps.7$z.fst) +0.7), fill = "red") +
#   geom_rect(data = genes.7, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(snps.7$z.fst) +0.7 , ymax= Inf, fill= type), color = "black" )+
#   scale_fill_manual(name = "", values = c("black", "white")) +
#   coord_cartesian( x = c(window.7$starts/1e6+ .075, window.7$ends/1e6 + .03)) +
#   labs(x = NULL, y =  "ZFst", subtitle = "CD209 antigen-like protein A") +
#   theme_classic() +
#   theme(legend.position = "none",
#         legend.background = element_rect(fill = "transparent"))

t.7 <- ggplot() +
  geom_rect(data= missense.gene.7, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.tajD.7, aes(x = BIN_START/1e6, TajimaD), color = "darkorange2", size =.75) +
  geom_line(data = LAO.tajD.7, aes(x = BIN_START/1e6, TajimaD), color = "blue4", size = .75) +
  coord_cartesian( x = c(window.7$starts/1e6, window.7$ends/1e6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL , y =  "Tajima's D") +
  theme_classic()

h.7 <- ggplot() +
  geom_rect(data= missense.gene.7, aes(xmin = start/1e6, xmax= end/1e6, ymin = -Inf, ymax = Inf, group = start), alpha = 0.25) +
  geom_line(data = STLO.het.7, aes(x =start/1e6, mean.hz, color = "GL Odd"), size =.75) +
  geom_line(data = LAO.het.7, aes(x = start/1e6, mean.hz, color = "BC Odd"), size =.75) +
  coord_cartesian( x = c(window.7$starts/1e6, window.7$ends/1e6), y = c(0,.7)) +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4", 
                                           "GL Odd" = "darkorange2"
  ))  +
  labs(x = paste("Chromosome", window.7$c, "position (Mbp)"), y =  "Heterozygosity") +
  theme_classic() +
  theme(legend.position = c(.1,.9),
        legend.background = element_rect(fill = "transparent"))

ow.7 <- s.7/t.7/h.7

### plot all genes (removed s.7 because it occurs in same window as s.1 (per2))

#ag <- (s.2/s.4/s.5/s.6/s.3) + plot_annotation(title = "c)")
ag <- (s.2/s.3) + plot_annotation(title = "c)")
#ag <- (s.2 | s.4 | s.5 | s.6 |s.3) + plot_annotation(title = "c)")
#ag <- ag & theme_classic(base_size = 8) + theme(legend.position = "none") 
ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/all.missense.genes.pdf", ag,  height = 5, width = 3.75, units = "in", dpi  = 600, )


# Extras ------------------------------------------------------------------




wind.adj <-50000

### plot gene w 4 missense

gene1 <- Zfst.missense[1,1]

gene1.pos <- Zfst.gene.pos[str_which(string = Zfst.gene.pos$attributes, pattern = as.character(gene1)),]

gene1.snps <- Zfst.snps[Zfst.snps$CHROM == gene1.pos[1,]$seqid & Zfst.snps$POS >= gene1.pos[1,]$start -wind.adj & Zfst.snps$POS <= gene1.pos[1,]$end +wind.adj,]         

p.g1 <- ggplot(gene1.snps) +
  geom_point(aes(x = POS/1e6, y = z.fst), size = 0.5) +
  geom_point(data = gene1.snps[gene1.snps$z.fst >=4,], aes(x = POS/1e6, y = z.fst), size = .5, color = "red") +
  geom_rect(data = gene1.pos, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(gene1.snps$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  labs(x ="Position (Mbp)", y = "ZFst", subtitle = "CD209 antigen-like protein A") +
  theme_classic()


### plot gene w 3 missense



gene2 <- Zfst.missense[2,1]

gene2.pos <- Zfst.gene.pos[str_which(string = Zfst.gene.pos$attributes, pattern = as.character(gene2)),]

gene2.snps <- Zfst.snps[Zfst.snps$CHROM == gene2.pos[1,]$seqid & Zfst.snps$POS >= gene2.pos[1,]$start -wind.adj & Zfst.snps$POS <= gene2.pos[1,]$end +wind.adj,]         

p.g2 <-ggplot(gene2.snps) +
  geom_point(aes(x = POS/1e6, y = z.fst), size = 0.5) +
  geom_point(data = gene2.snps[gene2.snps$z.fst >=5,], aes(x = POS/1e6, y = z.fst), size = .5, color = "red") +
  geom_rect(data = gene2.pos, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(gene2.snps$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  labs(x ="Position (Mbp)", y = "ZFst", subtitle = "period circadian protein homolog 2-like") +
  theme_classic()

### plot puberty genes

gene4 <- Zfst.missense[4,1]

gene4.pos <- Zfst.gene.pos[str_which(string = Zfst.gene.pos$attributes, pattern = as.character(gene4)),]

gene4.snps <- Zfst.snps[Zfst.snps$CHROM == gene4.pos[1,]$seqid & Zfst.snps$POS >= gene4.pos[1,]$start -wind.adj & Zfst.snps$POS <= gene4.pos[1,]$end +wind.adj,]         

p.g4 <-ggplot(gene4.snps) +
  geom_point(aes(x = POS/1e6, y = z.fst), size = 0.5) +
  geom_point(data = gene4.snps[gene4.snps$z.fst >=4,], aes(x = POS/1e6, y = z.fst), size = .5, color = "red") +
  geom_rect(data = gene4.pos, aes(xmin = start/1e6, xmax = end/1e6, ymin = max(gene4.snps$z.fst) +0.5 , ymax= Inf, fill= type), color = "black" )+
  scale_fill_manual(name = "", values = c("black", "white")) +
  labs(x ="Position (Mbp)", y = "ZFst", subtitle = "gnrhr4:\ngonadotropin releasing hormone receptor 4") +
  theme_classic()

p.g1/p.g2/p.g4
              

### genes near per2
Zfst.gene.pos.coding[str_which(Zfst.gene.pos.coding$attributes,pattern = "LOC124009266" ),]

mean(raw.fst[raw.fst$CHROM ==22 & raw.fst$POS >= 14693882 & raw.fst$POS <=14769178, "WEIR_AND_COCKERHAM_FST"], na.rm = TRUE)

Zfst.gene.pos.coding[str_which(Zfst.gene.pos.coding$attributes,pattern = "LOC124009839" ),]

mean(raw.fst[raw.fst$CHROM ==22 & raw.fst$POS >= 14759544 & raw.fst$POS <=14761621, "WEIR_AND_COCKERHAM_FST"], na.rm = TRUE)

Zfst.gene.pos.coding[str_which(Zfst.gene.pos.coding$attributes,pattern = "LOC124010088" ),]

mean(raw.fst[raw.fst$CHROM ==22 & raw.fst$POS >= 14763053 & raw.fst$POS <=14769178, "WEIR_AND_COCKERHAM_FST"], na.rm = TRUE)

Zfst.gene.pos.coding[str_which(Zfst.gene.pos.coding$attributes,pattern = "lpar5a" ),]

mean(raw.fst[raw.fst$CHROM ==22 & raw.fst$POS >= 14769294 & raw.fst$POS <=14770818, "WEIR_AND_COCKERHAM_FST"], na.rm = TRUE)


  
### function
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


### Count number of eigenGWAS SNPs 
egwas.snps <- read_tsv("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/eigenGWAS/LAO_STLO_samples.1.egwas")
egwas.snps$logpval <- -log10(egwas.snps$P)
egwas.snps <- egwas.snps[egwas.snps$logpval >= 6.64,]


eg.snps.app <- NULL
for (z in 1:nrow(outliers.winds.corr)){
  wind <- outliers.winds.corr[z,]
  w.c <- wind$c
  w.s <- wind$starts
  w.e <- wind$ends
  eg.snps <- nrow(na.omit(egwas.snps[egwas.snps$CHR == w.c & egwas.snps$BP >= w.s & egwas.snps$BP <= w.e, ]))
  
  egwas.window.counts <- cbind(wind, snps = eg.snps)
  eg.snps.app <- rbind( eg.snps.app, egwas.window.counts)
}

sum(eg.snps.app$snps)
