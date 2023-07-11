# Header ------------------------------------------------------------------
# This script calculates some basic stats and provides the plots for Fig. 2

setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon/")


library(tidyverse); library(ggdist)

# Heterozygosity estimates ------------------------------------------------

indv.het <- read_csv("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/indv.het.csv")

indv.het$Pop <- recode(indv.het$Pop, 
                       "LAO" = "BC Odd",
                       "LAE" = "BC Even",
                       "STLO" = "GL Odd",
                       "STLE" = "GL Even",
                       "STLO3" = "GL Odd 3")

het <- ggplot(indv.het, aes(y = Pop, x = het, fill =Pop, color = Pop)) +
  geom_point( size = 1.3, alpha = .25, position = position_jitter(width = .001)) +
  stat_halfeye( adjust = .5, width = .3, .width = 0, justification = -1, point_colour = NA, scale = 0.5, height = .75) +
  geom_boxplot(aes(middle =  mean(het)),width = .15,  alpha = .25,  outlier.shape = NA, ) +
  labs(y = "", x = "Heterozygosity", title = "b)") +
  scale_fill_manual(values = c("darkgreen","blue4",  "orange",  "darkorange2", "darkorange4")) +
  #scale_fill_manual(values = c("orange", "darkorange2", "darkorange4",  "darkgreen", "blue2" )) +
  scale_color_manual(values = c("darkgreen","blue4",  "orange",  "darkorange2", "darkorange4")) +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  lims(x = c(0.1,0.35)) +
  theme_classic() +
  theme(legend.position = "none")

#summary table by sample group
indv.het %>% 
  group_by(Pop) %>% 
  summarise(mean.het = mean(het), sd.het = sd(het)) %>% 
  arrange(desc(mean.het))

# find lowest individuals
indv.het %>%  dplyr::select(Pop, Sample, id, het)  %>% arrange(het)

# summary wihout those individuals
indv.het %>% 
  filter(het >= 0.12) %>% 
  group_by(Pop) %>% 
  summarise(mean.het = mean(het), sd.het = sd(het)) %>% 
  arrange(desc(mean.het))

# Windowed Heterozygosity -------------------------------------------------

### sccript to calulate windowed het from the big, joint vcf


# all.pops.vcf <- read.vcf("./data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/05_MAF/joint/GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.indv80.gntyp80.MAF005.recode.newnames.newchr.vcf", convert.chr = TRUE)
# 
# all.pops.vcf@ped$famid <-all.pops.vcf@ped %>% 
#   tibble() %>% 
#   separate(famid, sep = "_", into = c("pop", "ind")) %>% 
#   select("pop")
# 
# chr10.inv <- all.pops.vcf[,all.pops.vcf@snps$chr ==10 & all.pops.vcf@snps$pos >= 6.5e6 & all.pops.vcf@snps$pos <= 35e6]
# 
# ggplot() +
#   geom_histogram(data = chr10.inv@ped, aes(x = hz, color = famid$pop), fill = "white") +
#   facet_wrap(~famid$pop, ncol = 1) +
#   theme_classic()
# 
# # creat list of samples and then subset vcfs
# STLO.samples <- all.pops.vcf@ped$id[str_which(string = all.pops.vcf@ped$id, pattern = "STLO_")]
# STLO3.samples <- all.pops.vcf@ped$id[str_which(string = all.pops.vcf@ped$id, pattern = "STLO3_")]
# STLE.samples <- all.pops.vcf@ped$id[str_which(string = all.pops.vcf@ped$id, pattern = "STLE_")]
# LAO.samples <-  all.pops.vcf@ped$id[str_which(string = all.pops.vcf@ped$id, pattern = "LAO_")]
# LAE.samples <-  all.pops.vcf@ped$id[str_which(string = all.pops.vcf@ped$id, pattern = "LAE_")]
# 
# STLO.vcf <- select.inds(x = all.pops.vcf, condition = id %in% STLO.samples)
# STLO3.vcf <- select.inds(x = all.pops.vcf, condition = id %in% STLO3.samples)
# STLE.vcf <- select.inds(x = all.pops.vcf, condition = id %in% STLE.samples)
# LAO.vcf <- select.inds(x = all.pops.vcf, condition = id %in% LAO.samples)
# LAE.vcf <- select.inds(x = all.pops.vcf, condition = id %in% LAE.samples)


### function for windowed het across entire genome

wind.het <- function(vcf, step){
  chroms <- unique(vcf@snps$chr)
  
  OUT.all <- NULL
  for (c in 1:length(chroms)){
    
    chr <- chroms[c]
    # get snps slot and subset to chr
    data <- vcf
    snps <- data@snps
    snps <- snps[snps$chr == chr,]
    
    # calculate sliding windows
    
    window.start <- snps$pos[1] 
    window.end   <- snps$pos[length(snps[, 1])]
    
    wstarts <- seq(from = window.start, to = window.end, by = step)
    wends   <- wstarts + step
    
    #loop to calc over windows 
    OUT5 <- NULL
    for(w in 1:length(wstarts)){
      window   <- wstarts[w]:wends[w]
      snps2     <- snps[which(snps$pos >= wstarts[w] & snps$pos < wends[w]), ]
      hz <- ifelse(chr != 23, mean(snps2$hz), sum(snps2$N1)/sum(c(snps2$N0,snps2$N1,snps2$N2))) # ifelse for chr23 and calculate windowed hz
      mean.hz <- hz
      out      <- cbind(chr,w, wstarts[w], wends[w], mean.hz)
      OUT5     <- rbind(OUT5, out)
    }
    OUT.all <-rbind(OUT.all, OUT5)
    
  }
  
  windowed.hz <- data.frame(OUT.all)
  windowed.hz[,3] <- as.numeric(windowed.hz[,3])
  windowed.hz[,4] <- as.numeric(windowed.hz[,4])
  windowed.hz[,5] <- as.numeric(windowed.hz[,5])
  colnames(windowed.hz)[3:4] <- c("start", "end")
  
  return(windowed.hz)
  
}

# step.size <- 2.5e6
# 
# STLO.het <- wind.het(vcf = STLO.vcf, step = step.size)
# STLE.het <- wind.het(vcf = STLE.vcf, step = step.size)
# STLO3.het <- wind.het(vcf = STLO3.vcf, step = step.size)
# LAO.het <- wind.het(vcf = LAO.vcf, step = step.size)
# LAE.het <- wind.het(vcf = LAE.vcf, step = step.size)
# 
# 
# write.table(x = STLO.het, "./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/STLO.2.5MBwindowed.het.txt", sep = "\t")
# write.table(x = STLO3.het, "./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/STLO3.2.5MBwindowed.het.txt", sep = "\t")
# write.table(x = STLE.het, "./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/STLE.2.5MBwindowed.het.txt", sep = "\t")
# write.table(x = LAO.het, "./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/LAO.2.5MBwindowed.het.txt", sep = "\t")
# write.table(x = LAE.het, "./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/LAE.2.5MBwindowed.het.txt", sep = "\t")

STLO.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/STLO.2.5MBwindowed.het.txt")
STLO3.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/STLO3.2.5MBwindowed.het.txt")
STLE.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/STLE.2.5MBwindowed.het.txt")
LAO.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/LAO.2.5MBwindowed.het.txt")
LAE.het <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/LAE.2.5MBwindowed.het.txt")


het.g.raw <- ggplot() + 
  geom_line(data = LAO.het, aes(x = start, y = mean.hz , color = "BC Odd" ), size = 0.5) +
  geom_line(data = LAE.het, aes(x = start, y = mean.hz , color = "BC Even" ), size = 0.5) +
  geom_line(data = STLO.het, aes(x = start, y = mean.hz, color= "GL Odd"), size = 0.5) +
  geom_line(data = STLO3.het, aes(x = start, y = mean.hz , color = "GL Odd 3"), size = 0.5) +
  geom_line(data = STLE.het, aes(x = start, y = mean.hz , color = "GL Even" ), size = 0.5) +
  facet_grid(~factor(chr), switch = "x", scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "", values = c("BC Odd" = "Blue4",
                                           "BC Even" = "darkgreen",
                                           "GL Odd" = "darkorange2",
                                           "GL Even" = "orange",
                                           "GL Odd 3" = "darkorange4")) +
  labs(x = "Chromosome", y = "Heterozygosity") +
  theme_classic( base_size = 24) +
  theme(legend.position="top", 
        strip.background = element_rect(colour="white", fill="white", linewidth=1.5, linetype="solid"), #formatting for facets
        panel.spacing.x=unit(.1, "lines"),
        axis.title.x=element_blank(), # remove x axis 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Supplemental Materials/gwide.het.all.pdf", het.g.raw, height = 8, width = 23, units = "in", dpi = 600 )

## without extreme het values
# set up a box for inversion
inversion.box <- STLO.het[1,]
inversion.box$chr <- 10
inversion.box[1,"start"] <-6403689
inversion.box[1,"end"] <-35403689



het.g <- ggplot() + 
  geom_rect(data = inversion.box, aes(xmin = start, xmax = end, ymin = -Inf, ymax= Inf ), fill = "grey20", alpha = 0.2) +
  geom_line(data = STLO.het[STLO.het$mean.hz <=.35,], aes(x = start, y = mean.hz, color= "GL Odd"), linewidth = 0.5) +
  # geom_line(data = STLO3.het, aes(x = start, y = mean.hz) , color = "red", linewidth = 0.5) +
  # geom_line(data = STLE.het, aes(x = start, y = mean.hz) , color = "purple", linewidth = 0.5) +
  geom_line(data = LAO.het[ LAO.het$mean.hz <=.35,], aes(x = start, y = mean.hz, color = "BC Odd"), linewidth = 0.5) +
  #geom_line(data = LAE.het, aes(x = start, y = mean.hz) , color = "darkgreen", linewidth = 0.5) +
  scale_color_manual(values = c("BC Odd" = "blue4",
                                "GL Odd" = "darkorange"),
                     name = NULL) +
  facet_grid(~chr, switch = "x", scales = "free_x", space = "free_x") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Chromosome", y = "Heterozygosity", title = "a)") +
  theme_classic() +
   theme(legend.position="top", # no legened
        strip.background = element_rect(colour="white", fill="white", linewidth=1.5, linetype="solid"), #formatting for facets
        panel.spacing.x=unit(.1, "lines"),
        axis.title.x=element_blank(), # remove x axis 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/CH 3/Figures/Fig. 2a.pdf", het.g, height = 2, width = 8.5, units = "in", dpi = 600)

### summary het
inv.start <- 6403689
inv.end <- 35403689

chr10.het <-data.frame(matrix(data = NA, nrow = 5, ncol = 3))
colnames(chr10.het) <- c("pop", "non.ROI", "ROI")  
chr10.het$pop <- c("BC Odd", "BC Even", "GL Odd", "GL Odd 3", " GL Even")

chr10.het[1,"ROI"] <- mean(LAO.het[LAO.het$chr ==10 & LAO.het$start >= inv.start & LAO.het$end <= inv.end,"mean.hz"])
chr10.het[1,"non.ROI"] <-mean(c(LAO.het[LAO.het$chr ==10 & LAO.het$start <= inv.start,"mean.hz"], LAO.het[LAO.het$chr ==10 & LAO.het$end >= inv.end,"mean.hz"]))

chr10.het[2,"ROI"] <- mean(LAE.het[LAE.het$chr ==10 & LAE.het$start >= inv.start & LAE.het$end <= inv.end,"mean.hz"])
chr10.het[2,"non.ROI"] <-mean(c(LAE.het[LAE.het$chr ==10 & LAE.het$start <= inv.start,"mean.hz"], LAE.het[LAE.het$chr ==10 & LAE.het$end >= inv.end,"mean.hz"]))

chr10.het[3,"ROI"] <- mean(STLO.het[STLO.het$chr ==10 & STLO.het$start >= inv.start & STLO.het$end <= inv.end,"mean.hz"])
chr10.het[3,"non.ROI"] <- mean(c(STLO.het[STLO.het$chr ==10 & STLO.het$start <= inv.start,"mean.hz"], STLO.het[STLO.het$chr ==10 & STLO.het$end >= inv.end,"mean.hz"]))

chr10.het[4,"ROI"] <- mean(STLO3.het[STLO3.het$chr ==10 & STLO3.het$start >= inv.start & STLO3.het$end <= inv.end,"mean.hz"])
chr10.het[4,"non.ROI"] <- mean(c(STLO3.het[STLO3.het$chr ==10 & STLO3.het$start <= inv.start,"mean.hz"], STLO3.het[STLO3.het$chr ==10 & STLO3.het$end >= inv.end,"mean.hz"]))

chr10.het[5,"ROI"] <- mean(STLO3.het[STLO3.het$chr ==10 & STLO3.het$start >= inv.start & STLO3.het$end <= inv.end,"mean.hz"])
chr10.het[5,"non.ROI"] <- mean(c(STLO3.het[STLO3.het$chr ==10 & STLO3.het$start <= inv.start,"mean.hz"], STLO3.het[STLO3.het$chr ==10 & STLO3.het$end >= inv.end,"mean.hz"]))


chr10.het[,2:3] <- round(chr10.het[,2:3], 3)

chr10.het
# Total SNPs --------------------------------------------------------------

### these data are from nonMAF filtered population-specific VCFs with invariant sites removed (results of filter_invariantSNPs_byPop.sh file)

SNP.counts <- data.frame(matrix(NA, nrow = 5, ncol = 2))

colnames(SNP.counts) <- c("Population", "SNP.count")

SNP.counts$Population <- c("GL Odd", "GL Odd 3", "GL Even", "BC Odd", "BC Even")

SNP.counts$SNP.count <-c(5080474,5362161,4996010,7337647,7167282)

# SNP percent change (BC Odd to GL Odd)
(7337647-5080474)/7337647

((SNP.counts[4,2]-SNP.counts[3,2])/SNP.counts[4,2])*100
((SNP.counts[4,2]-SNP.counts[2,2])/SNP.counts[4,2])*100
((SNP.counts[4,2]-SNP.counts[1,2])/SNP.counts[4,2])*100


SNPS <- ggplot(SNP.counts, aes(x = Population, y = SNP.count/1e6, fill = Population)) +
  geom_col() +
  scale_fill_manual(values = c("darkgreen","blue4",  "orange",  "darkorange2", "darkorange4")) +
  labs(x = NULL, y = "Number of SNPs\n(millions)", title = "c)") +
  theme_classic() +
  lims(y = c(0,8)) +
  theme(legend.position = "none")

# SNPs per window ---------------------------------------------------------

windowed.snps <- read_tsv("~/Downloads/snps_per_chromosome.txt")

colnames(windowed.snps)[1] <-  "Pop"

windowed.snps$Pop <- recode(windowed.snps$Pop, "1" = "BC Odd", "2" = "GL Odd")

windowed.snps$chromosome <- recode(windowed.snps$chromosome,
                                   "NC_060173.1" = 1 ,
                                   "NC_060174.1" = 2 ,
                                   "NC_060175.1" = 3 ,
                                   "NC_060176.1" = 4 ,
                                   "NC_060177.1" = 5 ,
                                   "NC_060178.1" = 6 ,
                                   "NC_060179.1" = 7 ,
                                   "NC_060180.1" = 8 ,
                                   "NC_060181.1" = 9 ,
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
                                   "NC_060198.1" = 26)


windowed.snps %>% group_by(chromosome, Pop) %>%
  summarise(avg_snps = mean(total_n_snps)) %>% 
  ggplot(aes(x = chromosome, y = avg_snps, fill = Pop)) +
  geom_col(position="dodge") +
  scale_fill_manual(name = "", values = c(
    "GL Odd" = "darkorange2",
    "BC Odd" = "blue4")) + 
  labs(x = "Chromosome", y = "Total SNPs per 100 Kbp", subtitle = "37% loss of SNPs") +
  theme_classic()

windowed.snps %>% group_by(Pop) %>% 
  summarise(avg_snps = mean(total_n_snps))



# function to resize window snps (here 25 = 25 windows of 100 kbp or 2.mbp)
resize.window.snps <- function(window){
  OUT <- NULL
 for(p in unique(windowed.snps$Pop)){
   
   sub.df <- windowed.snps[windowed.snps$Pop == p,]

 OUT.1 = NULL
  for(c in unique(sub.df$chromosome)){
    sub.c.df <- sub.df[sub.df$chromosome == c,]
    
    # calculate sliding windows
    
    window.start <- 1
    window.end   <- nrow(sub.c.df)
    
    wstarts <- seq(from = window.start, to = window.end, by = window)
    wends   <- wstarts + window
  OUT.2 <- NULL
  for(w in 1:length(wstarts)){
    sub.c.w <- sub.c.df[wstarts[w]:wends[w],]
    mean_n_snps <- sum(sub.c.w$total_n_snps)
    temp.row <- cbind(sub.c.w[1,1:3], mean_n_snps = mean_n_snps)
    OUT.2<- rbind(OUT.2, temp.row )
    }
  OUT.1 <- rbind(OUT.1, OUT.2)
  }
 OUT <- rbind(OUT, OUT.1)
 
 }
  OUT <- data.frame(OUT)
  return(OUT)
} 
  
windowed.snps.2.5mbp <- resize.window.snps(25)
  

w.snps <- ggplot(windowed.snps.2.5mbp) +
  #geom_point(aes(x = window_start_position, y = mean_n_snps, group = Pop, color = Pop), size = 1) +
  geom_line(aes(x = window_start_position, y = mean_n_snps, group = Pop, color = Pop)) +
  facet_grid(~factor(chromosome), switch = "x", scales = "free_x", space = "free_x") +
  scale_color_manual(name = "", values = c(
    "GL Odd" = "darkorange2",
    "BC Odd" = "blue4")) + 
  labs(x = "Chromosome", y = "Number of SNPs", title = "a)") +
  lims(y = c(0,17500)) +
  theme_classic() +
  theme(legend.position="none", 
        strip.background = element_rect(colour="white", fill="white", linewidth=1.5, linetype="solid"), #formatting for facets
        panel.spacing.x=unit(.1, "lines"),
        axis.title.x=element_blank(), # remove x axis 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

windowed.snps.2.5mbp %>% 
  group_by(Pop) %>% 
  summarise(avg_SNPs = mean(mean_n_snps, na.rm = TRUE))

# Pairwise Fst ------------------------------------------------------------
pairwise.fst <-read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/SNPrelate/pairwise.fst.txt", header = TRUE)

pairwise.fst$pop1 <- recode(pairwise.fst$pop1, 
                            "LAO" = "BC Odd",
                            "LAE" = "BC Even",
                            "STLO" = "GL Odd",
                            "STLE" = "GL Even",
                            "STLO3" = "GL Odd 3")

pairwise.fst$pop2 <- recode(pairwise.fst$pop2, 
                            "LAO" = "BC Odd",
                            "LAE" = "BC Even",
                            "STLO" = "GL Odd",
                            "STLE" = "GL Even",
                            "STLO3" = "GL Odd 3")

p.fst <- ggplot(pairwise.fst) +
  geom_tile(aes(x= pop1, y = pop2, fill =  fst), color = "white") +
  scale_fill_gradient(limits = c(0.00007, .17), low = "blue4", high = "red", name = bquote(paste(italic('F'['ST']*' ')))) + 
  geom_text(aes(x= pop1, y = pop2, label = fst), color = "white", size = 3) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(.9,.3))

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/CH 3/Figures/Fig. 1c.pdf", p.fst, height = 4, width = 4, units = "in", dpi = 600)

  

# Combine plots -----------------------------------------------------------

fig2 <- w.snps /  (het + SNPS ) +plot_annotation(title = 'Figure 2')

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/Fig 2.pdf", fig2,  height = 4.5, width = 8.5, units = "in", dpi = 600)


# Extra Plots -------------------------------------------------------------


#### TajD 

# LAO_tajD <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/TajD/noMAF/", header = TRUE, sep = "\t")
# 
# LAE_tajD <- read.table("./pop_gen/TajD/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD5MB.Tajima.D", header = TRUE, sep = "\t")
# 
# STLO_tajD <- read.table("./pop_gen/TajD/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD5MB.Tajima.D", header = TRUE, sep = "\t")
# 
# STLO3_tajD <- read.table("./pop_gen/TajD/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD5MB.Tajima.D", header = TRUE, sep = "\t")
# 
# STLE_tajD <- read.table("./pop_gen/TajD/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD5MB.Tajima.D", header = TRUE, sep = "\t")
# 
# TajD <- ggplot() + 
#   geom_line(data = STLO_tajD, aes(x = BIN_START, y = TajimaD), color= "darkorange2") +
#   #geom_line(data = STLO3_tajD, aes(x = BIN_START, y = TajimaD) , color = "darkorange4") +
#   geom_line(data = STLE_tajD, aes(x = BIN_START, y = TajimaD) , color = "purple") +
#   #geom_line(data = LAO_tajD, aes(x = BIN_START, y = TajimaD) , color = "blue4") +
#   #geom_line(data = LAE_tajD, aes(x = BIN_START, y = TajimaD) , color = "darkgreen") +
#   facet_grid(~as.factor(CHROM), switch = "x",scales = "free_x", space = "free_x") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Chromosome", y = "Tajima's D") +
#   theme_classic() +
#   theme(legend.position="none", # no legened
#         strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"), #formatting for facets
#         panel.spacing.x=unit(.1, "lines"),
#         axis.title.x=element_blank(), # remove x axis 
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())

# STLO_tajD %>% 
#   arrange(TajimaD) %>% 
#   head
# 
# 
# mean(STLO_tajD$TajimaD)
# min(STLO_tajD$TajimaD)
# mean(LAO_tajD$TajimaD)
# mean(LAE_tajD$TajimaD)




#### Pi 
# LAO_pi <- read.table("./pop_gen/pi/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB.windowed.pi", header = TRUE, sep = "\t")
# 
# LAE_pi <- read.table("./pop_gen/pi/LAE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB.windowed.pi", header = TRUE, sep = "\t")
# 
# STLO_pi <- read.table("./pop_gen/pi/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB.windowed.pi", header = TRUE, sep = "\t")
# 
# STLO3_pi <- read.table("./pop_gen/pi/STLO3_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB.windowed.pi", header = TRUE, sep = "\t")
# 
# STLE_pi <- read.table("./pop_gen/pi/STLE_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi5MB.windowed.pi", header = TRUE, sep = "\t")
# 
# pi <- 
#   
#   
#   STLO_pi %>% 
#   filter(CHROM == 10) %>% 
#   arrange(PI) 
# 
# STLO_pi %>% 
#   filter(CHROM == 17) %>% 
#   arrange(PI) 


# ### Chrom 10 TajD and PI
# LAO_pi1mb <- read.table("./pop_gen/pi/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi1MB.windowed.pi", header = TRUE, sep = "\t")
# 
# STLO_pi1mb <- read.table("./pop_gen/pi/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.pi1MB.windowed.pi", header = TRUE, sep = "\t")
# 
# LAO_tajD1mb <- read.table("./pop_gen/TajD/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD1MB.Tajima.D", header = TRUE, sep = "\t")
# 
# STLO_tajD1mb <- read.table("./pop_gen/TajD/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD1MB.Tajima.D", header = TRUE, sep = "\t")
# 
# 
# 
# z
# pi_10 <- ggplot() + 
#   geom_line(data = STLO_pi1mb[STLO_pi1mb$CHROM ==10,], aes(x = BIN_START, y = PI), color= "darkorange2") +
#   geom_line(data = LAO_pi1mb[LAO_pi1mb$CHROM == 10,], aes(x = BIN_START, y = PI) , color = "blue4") +
#   facet_grid(~as.factor(CHROM), switch = "x",scales = "free_x") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Chromosome", y = "Pi") +
#   theme_classic() +
#   theme(legend.position="none", # no legened
#         strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"), #formatting for facets
#         panel.spacing.x=unit(.1, "lines"),
#         
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# tajdD_10 <- ggplot() + 
#   geom_line(data = STLO_tajD1mb[STLO_tajD1mb$CHROM ==10,], aes(x = BIN_START, y = TajimaD), color= "darkorange2") +
#   geom_line(data = LAO_tajD1mb[LAO_tajD1mb$CHROM == 10,], aes(x = BIN_START, y = TajimaD) , color = "blue4") +
#   facet_grid(~as.factor(CHROM), switch = "x",scales = "free_x") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Chromosome", y = "Tajima's D") +
#   theme_classic() +
#   theme(legend.position="none", # no legened
#         strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"), #formatting for facets
#         panel.spacing.x=unit(.1, "lines"),
#         
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# pi_10 / tajdD_10
# 
# 
# 
