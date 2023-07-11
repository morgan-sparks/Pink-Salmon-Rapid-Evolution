# Header ------------------------------------------------------------------

setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon/data/")
setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/archive/from_fortress/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/08_genotypeGVCFs/")

### libraries 

library(tidyverse); library(patchwork)


# Data --------------------------------------------------------------------

zfst.snps <- read.table("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAOvsSTLO.Zfst.snps.txt", header = TRUE)

zfst.windows <- read.table("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers/LAOvSTLOoutlier_windows_all_individuals.txt", header = TRUE)

zfst.100kb <- read.table("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAO_STLO.zfst.windows.txt", header = TRUE)

zfst.windows.corrected <- read.table("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers/LAOvSTLO_outlierwindows_corrected.txt", header = TRUE, sep = "\t")

egwas.snps <- read_tsv("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/eigenGWAS/LAO_STLO_samples.1.egwas")

egwas.snps$logpval <- -log10(egwas.snps$P)

STLO.tajd.10kb <- read.table("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/TajD/STLO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB.Tajima.D", header = TRUE)
LAO.tajd.10kb <- read.table("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/TajD/LAO_hardfiltered_snps.filt.biallelic.indiv80.MMCN.MAF5.renamed.TajimaD10KB.Tajima.D", header = TRUE)

STLO.het <-read.table("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/STLO.10KBwindowed.het.txt", header =  TRUE)
LAO.het <-read.table("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/heterozygosity/LAO.10KBwindowed.het.txt", header =  TRUE)

pcadapt.snps <- df
# Plot windows ------------------------------------------------------------

pdf("~/Downloads/outliers.pdf")
for (w in 1:nrow(zfst.windows.corrected)){
  single.window <- zfst.windows.corrected[w,]
  chr <- single.window$c #window chrom
  w.s <- single.window$starts # window start
  w.e <- single.window$ends # window end
  
  ### snps
  window.snps <- zfst.snps[zfst.snps$CHROM == chr & 
                             zfst.snps$POS >= w.s  & 
                             zfst.snps$POS <= w.e & 
                             zfst.snps$z.fst >= 0,]
  window.snps.4zfst <- window.snps[window.snps$z.fst >= 4,]
  
  # window snps +/- adjusted interval
  
  w.adj <- .5e6
  window.snps.wide <- zfst.snps[zfst.snps$CHROM == chr & 
                                  zfst.snps$POS >= w.s - w.adj  & 
                                  zfst.snps$POS <= w.e +w.adj & 
                                  zfst.snps$z.fst >= 0,]
  

  # window snps +/- adjusted interval
  
  window.snps.wide <-  zfst.snps[ zfst.snps$CHROM == chr & 
                                    zfst.snps$POS >= w.s - w.adj  & 
                                    zfst.snps$POS <= w.e +w.adj & 
                                    zfst.snps$z.fst >= 0,]
  
  ## egwas snps
  
  window.egwas <- egwas.snps[egwas.snps$CHR == chr & 
                               egwas.snps$BP >= w.s  & 
                               egwas.snps$BP <= w.e, ]
  
  window.egwas.logp8 <-window.egwas[window.egwas$logpval >= -log10(0.05/nrow(zfst.snps)), ]
  # wide egwas.snps
  
  window.egwas.wide <- egwas.snps[egwas.snps$CHR == chr & 
                                    egwas.snps$BP >= w.s - w.adj & 
                                    egwas.snps$BP <= w.e + w.adj,]
  
 
  # summary stats for interval
  
  STLO.tajd.window.wide <- STLO.tajd.10kb[STLO.tajd.10kb$CHROM == chr & 
                                            STLO.tajd.10kb$BIN_START >= w.s - w.adj & 
                                            STLO.tajd.10kb$BIN_START <= w.e + w.adj, ]
  
  LAO.tajd.window.wide <- LAO.tajd.10kb[LAO.tajd.10kb$CHROM == chr & 
                                            LAO.tajd.10kb$BIN_START >= w.s - w.adj & 
                                            LAO.tajd.10kb$BIN_START <= w.e + w.adj, ]
  
  STLO.window.het.wide <- STLO.het[STLO.het$chr == chr & 
                                            STLO.het$start >= w.s - w.adj & 
                                            STLO.het$start <= w.e + w.adj, ]
  
  LAO.window.het.wide <- LAO.het[LAO.het$chr == chr & 
                                          LAO.het$start >= w.s - w.adj & 
                                          LAO.het$start <= w.e + w.adj, ]
  
  
 p1 <- ggplot() +
    geom_point(data = window.snps.wide, aes(x = POS/1e6, y = z.fst), size = 0.5, color = "grey") +
    geom_point(data = window.snps, aes(x = POS/1e6, y = z.fst), size = 0.5, color = "black") +
    geom_point(data = window.snps.4zfst, aes(x = POS/1e6, y = z.fst), size = 0.5, color = "red") +
    geom_hline(yintercept = 4, linetype = "dashed") +
    labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "ZFst", subtitle = paste("Outlier window", w, "(length = ", w.e - w.s, "bp)")) +
    theme_classic()
  
 p2 <- ggplot() +
    geom_point(data = window.egwas.wide, aes(x = BP/1e6, y = logpval), size = 0.5, color = "grey") +
    geom_point(data = window.egwas, aes(x = BP/1e6, y = logpval), size = 0.5, color = "black") +
    geom_point(data = window.egwas.logp8, aes(x = BP/1e6, y = logpval), size = 0.5, color = "red") +
    geom_hline(yintercept = -log10(0.05/nrow(zfst.snps)), linetype = "dashed") +
    labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "-log10(P)") +
    theme_classic()
  
 p3 <-ggplot() +
    geom_line(data = STLO.tajd.window.wide, aes(x =BIN_START/1e6, y = TajimaD), color = "darkorange2") +
    geom_line(data = LAO.tajd.window.wide, aes(x =BIN_START/1e6, y = TajimaD), color = "blue4") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "Tajima's D\n(10 Kbp)") +
    theme_classic()
 
 p4 <-ggplot() +
   geom_line(data =  STLO.window.het.wide , aes(x =start/1e6, y = mean.hz, color = "GL Odd")) +
   geom_line(data =  LAO.window.het.wide , aes(x =start/1e6, y = mean.hz, color = "BC Odd"))+
   geom_hline(yintercept = 0, linetype = "dashed") +
   scale_color_manual( name = "", values = c("BC Odd" = "Blue4", 
                          "GL Odd" = "darkorange2")) +
   labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "Heterozygosity\n(10 Kbp)") +
   theme_classic() +
   theme(legend.position = "bottom")
  

#plot(p1/p2/p3/p4)

p <- p1/p2/p3/p4

assign(paste0("plot_", w), p)

plots[[w]] <- paste0("plot_", w)

   
}
#dev.off()


# Manhattan plot ----------------------------------------------------------

### make windows 
outliers.winds.corr <- read.table("./seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers/LAOvSTLO_outlierwindows_corrected.txt", header = TRUE)

misssense.windows <- outliers.winds.corr[c(3,32),]
misssense.windows$starts <- misssense.windows$starts- 5e6
misssense.windows$ends <- misssense.windows$ends + 5e6

per2.window <- outliers.winds.corr[c(31),]

per2.window$starts <- per2.window$starts - 5e6
per2.window$ends <-per2.window$ends + 5e6
### Zfst manhattan

chromosomes <- as.factor(unique(zfst.100kb$c))
chrom.cols <- c(rep(c("#000000" , "DarkViolet"), length(chromosomes)/2))


fst <- ggplot() +
  geom_rect(data = misssense.windows, aes(xmin = starts, xmax = ends, ymin = -Inf, ymax = Inf), fill = "grey", alpha = .5) +
  geom_rect(data = per2.window, aes(xmin = starts, xmax = ends, ymin = -Inf, ymax = Inf), fill = "dodgerblue", alpha = .5) +
  geom_point(data = zfst.100kb, aes(x = real.starts, y =z.fst, c, size = 0.1, shape = 16))+
  geom_point(data =  zfst.100kb[zfst.100kb$z.fst> 5, ],aes(x = real.starts, y =z.fst), size = 0.5, shape = 16, color = "red")+
  facet_grid(.~c, switch = "x", scales = "free_x", space = "free_x") +
  scale_y_continuous(limits=c(-1,9)) +
  scale_color_manual(values=chrom.cols) + 
  xlab(NULL)+
  ylab(bquote(paste('Z',italic('F'['ST']*' '), sep = ""))) +
  labs(subtitle = "a)", title = "BC Odd vs. GL Odd") +
  geom_hline(yintercept = 5, linetype = "dashed", size = 0.5) +
  theme_classic(base_size = 12) +
  theme(legend.position="none", # no legened
        strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"), #formatting for facets
        panel.spacing.x=unit(.1, "lines"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



### egwas manhattan

window.length <- 100000 #100kb
window.step   <- 50000  #50kb
chroms <- unique(egwas.snps$CHR)

egwas.snps.quick <- egwas.snps[, c("CHR", "BP", "logpval")]

OUT2 <- NULL
for(c in 1:length(chroms)){
  
  fst2 <- egwas.snps.quick[which(egwas.snps.quick$CHR == chroms[c]), ]
  
  # calculate z score
  
  window.start <- fst2$BP[1]
  window.end   <- fst2$BP[nrow(fst2)]
  
  wstarts <- seq(from = window.start, to = window.end, by = window.step)
  wends   <- wstarts + window.length
  
  OUT3 <- NULL
  for(w in 1:length(wstarts)){
    window   <- wstarts[w]:wends[w]
    fst3     <- fst2[which(fst2$BP >= wstarts[w] & fst2$BP < wends[w]), ]
    fsts     <- fst3$logpval
    mean.logpval <- mean(fsts, na.rm = TRUE)
    out      <- cbind(c, w, wstarts[w], wends[w], mean.logpval)
    OUT3     <- rbind(OUT3, out)
  }

  OUT2    <- rbind(OUT2, OUT3)
  
}

egwas.100kb.winds <-data.frame(OUT2)
colnames(egwas.100kb.winds)[3:4] <- c("start", "end")


chromosomes <- as.factor(unique(egwas.100kb.winds$c))
chrom.cols <- c(rep(c("#000000" , "DarkViolet"), length(chromosomes)/2))


cutoff <- -log10(0.01/nrow(egwas.100kb.winds))


man <- ggplot() +
  geom_rect(data = misssense.windows, aes(xmin = starts, xmax = ends, ymin = -Inf, ymax = Inf), fill = "grey", alpha = .5) +
  geom_rect(data = per2.window, aes(xmin = starts, xmax = ends, ymin = -Inf, ymax = Inf), fill = "dodgerblue", alpha = .5) +
  geom_point(data = egwas.100kb.winds, aes(x = start, y =mean.logpval, colour = as.factor(c)), size = 0.1, shape = 16)+
  geom_point(data =  egwas.100kb.winds %>% filter(mean.logpval > cutoff),aes(x = start, y =mean.logpval),  size = 0.5, shape = 16, color = "red")+
  facet_grid(.~c, switch = "x", scales = "free_x", space = "free_x") +
  scale_color_manual(values=chrom.cols) + 
  xlab("Chromosome")+
  ylab("-log10(P)") +
  labs(subtitle = "b)") +
  geom_hline(yintercept = cutoff, linetype = "dashed", size = 0.5) +
  theme_classic(base_size = 12) +
  theme(legend.position="none", # no legened
        strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"), #formatting for facets
        panel.spacing.x=unit(.1, "lines"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


fig.4.top <- fst/man


ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/fig.4.top.png", fig.4.top, width = 8.5, height = 3.5, dpi = 600)



# Plot Outliers for SI ----------------------------------------------------
# plot function to later use lapply
plot.outliers.SI <- function(w){
  
  single.window <- zfst.windows.corrected[w,]
  chr <- single.window$c #window chrom
  w.s <- single.window$starts # window start
  w.e <- single.window$ends # window end
  
  ### snps
  window.snps <- zfst.snps[zfst.snps$CHROM == chr & 
                             zfst.snps$POS >= w.s  & 
                             zfst.snps$POS <= w.e & 
                             zfst.snps$z.fst >= 0,]
  window.snps.4zfst <- window.snps[window.snps$z.fst >= 4,]
  
  # window snps +/- adjusted interval
  
  w.adj <- .5e6
  window.snps.wide <- zfst.snps[zfst.snps$CHROM == chr & 
                                  zfst.snps$POS >= w.s - w.adj  & 
                                  zfst.snps$POS <= w.e +w.adj & 
                                  zfst.snps$z.fst >= 0,]
  
  
  # window snps +/- adjusted interval
  
  window.snps.wide <-  zfst.snps[ zfst.snps$CHROM == chr & 
                                    zfst.snps$POS >= w.s - w.adj  & 
                                    zfst.snps$POS <= w.e +w.adj & 
                                    zfst.snps$z.fst >= 0,]
  
  ## egwas snps
  
  # window.egwas <- egwas.snps[egwas.snps$CHR == chr & 
  #                              egwas.snps$BP >= w.s  & 
  #                              egwas.snps$BP <= w.e, ]
  # 
  # window.egwas.logp8 <-window.egwas[window.egwas$logpval >= -log10(0.05/nrow(zfst.snps)), ]
  # # wide egwas.snps
  # 
  # window.egwas.wide <- egwas.snps[egwas.snps$CHR == chr & 
  #                                   egwas.snps$BP >= w.s - w.adj & 
  #                                   egwas.snps$BP <= w.e + w.adj,]
  # 
  ## pcadapt snps
  
  window.pcadapt <- pcadapt.snps[pcadapt.snps$chroms == chr & 
                               pcadapt.snps$positions >= w.s  & 
                               pcadapt.snps$positions <= w.e, ]
  
  window.pcadapt.logp3<-window.pcadapt[window.pcadapt$logqval >= 3,]

  # wide pcadapt.snps
  
  window.pcadapt.wide <- pcadapt.snps[pcadapt.snps$chroms == chr & 
                                    pcadapt.snps$positions >= w.s - w.adj & 
                                    pcadapt.snps$positions <= w.e + w.adj,]
  
  # summary stats for interval
  
  STLO.tajd.window.wide <- STLO.tajd.10kb[STLO.tajd.10kb$CHROM == chr & 
                                            STLO.tajd.10kb$BIN_START >= w.s - w.adj & 
                                            STLO.tajd.10kb$BIN_START <= w.e + w.adj, ]
  
  LAO.tajd.window.wide <- LAO.tajd.10kb[LAO.tajd.10kb$CHROM == chr & 
                                          LAO.tajd.10kb$BIN_START >= w.s - w.adj & 
                                          LAO.tajd.10kb$BIN_START <= w.e + w.adj, ]
  
  STLO.window.het.wide <- STLO.het[STLO.het$chr == chr & 
                                     STLO.het$start >= w.s - w.adj & 
                                     STLO.het$start <= w.e + w.adj, ]
  
  LAO.window.het.wide <- LAO.het[LAO.het$chr == chr & 
                                   LAO.het$start >= w.s - w.adj & 
                                   LAO.het$start <= w.e + w.adj, ]
  
  
  p1 <- ggplot() +
    geom_point(data = window.snps.wide, aes(x = POS/1e6, y = z.fst), size = 0.5, color = "grey") +
    geom_point(data = window.snps, aes(x = POS/1e6, y = z.fst), size = 0.5, color = "black") +
    geom_point(data = window.snps.4zfst, aes(x = POS/1e6, y = z.fst), size = 0.5, color = "red") +
    geom_hline(yintercept = 4, linetype = "dashed") +
    labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "ZFst", subtitle = paste("Outlier window", w, "(length = ", w.e - w.s, "bp)")) +
    theme_classic(base_size = 6)
  
  # p2 <- ggplot() +
  #   geom_point(data = window.egwas.wide, aes(x = BP/1e6, y = logpval), size = 0.5, color = "grey") +
  #   geom_point(data = window.egwas, aes(x = BP/1e6, y = logpval), size = 0.5, color = "black") +
  #   geom_point(data = window.egwas.logp8, aes(x = BP/1e6, y = logpval), size = 0.5, color = "red") +
  #   geom_hline(yintercept = -log10(0.05/nrow(zfst.snps)), linetype = "dashed") +
  #   labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "-log10(P)") +
  #   theme_classic(base_size = 6)
  
  p2 <- ggplot() +
    geom_point(data = window.pcadapt.wide, aes(x = positions/1e6, y = logqval), size = 0.5, color = "grey") +
    geom_point(data = window.pcadapt, aes(x = positions/1e6, y = logqval), size = 0.5, color = "black") +
    geom_point(data = window.pcadapt.logp3, aes(x = positions/1e6, y = logqval), size = 0.5, color = "red") +
    geom_hline(yintercept = 3, linetype = "dashed") +
    labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "-log10(Q)") +
    ylim(0, max(window.pcadapt.wide$logqval)) +
    theme_classic(base_size = 6)
  
  p3 <-ggplot() +
    annotate("rect", fill = "grey60", alpha = .25, 
             xmin = (w.s)/1e6, xmax = (w.e )/1e6,
             ymin = -Inf, ymax = Inf) +
    geom_line(data = STLO.tajd.window.wide, aes(x =BIN_START/1e6, y = TajimaD), color = "darkorange2") +
    geom_line(data = LAO.tajd.window.wide, aes(x =BIN_START/1e6, y = TajimaD), color = "blue4") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "Tajima's D\n(10 Kbp)") +
    theme_classic(base_size = 6)
  
  p4 <-ggplot() +
    annotate("rect", fill = "grey60", alpha = .25, 
             xmin = (w.s)/1e6, xmax = (w.e )/1e6,
             ymin = -Inf, ymax = Inf) +
    geom_line(data =  STLO.window.het.wide , aes(x =start/1e6, y = mean.hz, color = "GL Odd")) +
    geom_line(data =  LAO.window.het.wide , aes(x =start/1e6, y = mean.hz, color = "BC Odd"))+
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual( name = "", values = c("BC Odd" = "Blue4", 
                                              "GL Odd" = "darkorange2")) +
    labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "Heterozygosity\n(10 Kbp)") +
    theme_classic(base_size = 6) +
    theme(legend.position = "bottom")
  
  
  all.p <- p1/p2/p3/p4

  all.p
  
}

### lapply into a list of plots
myplots <- lapply(c(1:nrow(zfst.windows.corrected)), plot.outliers.SI)

### loop to plot 4 at a time into a pdf

i =1
pdf("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/SI_outliers.pdf", height = 11, width = 8.5)
pdf("~/Downloads/SI_outliers.pdf", height = 11, width = 8.5)
while( i <= nrow(zfst.windows.corrected)) {
  
  if(i < 30){
  myplot <- patchwork::wrap_plots(myplots[c(i, i+ 1, i + 2, i+3, i+4, i+5)], ncol=2, nrow = 3)
  plot(myplot)
  i = i + 6
  }else{
    myplot <- patchwork::wrap_plots(myplots[c(31,32,33,34)], ncol=2, nrow = 3)
    plot(myplot)
    i = i + 4
  }
  }
dev.off()


# Plot casein kinase I (phosphorylator for Per2) --------------------------

# middle of CK1 gene
CK1.length <-59121991 - 59086514
CK1.middle <- 59121991 -CK1.length/2

CK1.window <- data.frame(cbind(c = 13, starts = CK1.middle - 5e4, ends = CK1.middle + 5e4))



### Plot casein kinase I

single.window <- CK1.window
chr <- single.window$c #window chrom
w.s <- single.window$starts # window start
w.e <- single.window$ends # window end

### snps
window.snps <- zfst.snps[zfst.snps$CHROM == chr & 
                           zfst.snps$POS >= w.s  & 
                           zfst.snps$POS <= w.e & 
                           zfst.snps$z.fst >= 0,]
window.snps.4zfst <- window.snps[window.snps$z.fst >= 4,]

# window snps +/- adjusted interval

w.adj <- .5e6
window.snps.wide <- zfst.snps[zfst.snps$CHROM == chr & 
                                zfst.snps$POS >= w.s - w.adj  & 
                                zfst.snps$POS <= w.e +w.adj & 
                                zfst.snps$z.fst >= 0,]


# window snps +/- adjusted interval

window.snps.wide <-  zfst.snps[ zfst.snps$CHROM == chr & 
                                  zfst.snps$POS >= w.s - w.adj  & 
                                  zfst.snps$POS <= w.e +w.adj & 
                                  zfst.snps$z.fst >= 0,]

## egwas snps

window.egwas <- egwas.snps[egwas.snps$CHR == chr & 
                             egwas.snps$BP >= w.s  & 
                             egwas.snps$BP <= w.e, ]

window.egwas.logp8 <-window.egwas[window.egwas$logpval >= -log10(0.05/nrow(zfst.snps)), ]
# wide egwas.snps

window.egwas.wide <- egwas.snps[egwas.snps$CHR == chr & 
                                  egwas.snps$BP >= w.s - w.adj & 
                                  egwas.snps$BP <= w.e + w.adj,]


# summary stats for interval

STLO.tajd.window.wide <- STLO.tajd.10kb[STLO.tajd.10kb$CHROM == chr & 
                                          STLO.tajd.10kb$BIN_START >= w.s - w.adj & 
                                          STLO.tajd.10kb$BIN_START <= w.e + w.adj, ]

LAO.tajd.window.wide <- LAO.tajd.10kb[LAO.tajd.10kb$CHROM == chr & 
                                        LAO.tajd.10kb$BIN_START >= w.s - w.adj & 
                                        LAO.tajd.10kb$BIN_START <= w.e + w.adj, ]

STLO.window.het.wide <- STLO.het[STLO.het$chr == chr & 
                                   STLO.het$start >= w.s - w.adj & 
                                   STLO.het$start <= w.e + w.adj, ]

LAO.window.het.wide <- LAO.het[LAO.het$chr == chr & 
                                 LAO.het$start >= w.s - w.adj & 
                                 LAO.het$start <= w.e + w.adj, ]

CK1.start <- 59086514/1e6
CK1.end <- 59121991/1e6


p1 <- ggplot() +
  geom_rect(aes(xmin = CK1.start, xmax = CK1.end, ymin = -Inf, ymax = Inf), fill = "forestgreen", alpha = .5) +
  geom_point(data = window.snps.wide, aes(x = POS/1e6, y = z.fst), size = 0.5, color = "grey") +
  geom_point(data = window.snps, aes(x = POS/1e6, y = z.fst), size = 0.5, color = "black") +
  geom_point(data = window.snps.4zfst, aes(x = POS/1e6, y = z.fst), size = 0.5, color = "red") +
  geom_hline(yintercept = 4, linetype = "dashed") +
  labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "ZFst", title = "Casein Kinase 1 (dark green)") +
  theme_classic(base_size = 12)

p2 <- ggplot() +
  geom_rect(aes(xmin = CK1.start, xmax = CK1.end, ymin = -Inf, ymax = Inf), fill = "forestgreen", alpha = .5) +
  geom_point(data = window.egwas.wide, aes(x = BP/1e6, y = logpval), size = 0.5, color = "grey") +
  geom_point(data = window.egwas, aes(x = BP/1e6, y = logpval), size = 0.5, color = "black") +
  geom_point(data = window.egwas.logp8, aes(x = BP/1e6, y = logpval), size = 0.5, color = "red") +
  geom_hline(yintercept = -log10(0.05/nrow(zfst.snps)), linetype = "dashed") +
  labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "-log10(P)") +
  theme_classic(base_size = 12)

p3 <-ggplot() +
  geom_rect(aes(xmin = CK1.start, xmax = CK1.end, ymin = -Inf, ymax = Inf), fill = "forestgreen", alpha = .5) +
  geom_line(data = STLO.tajd.window.wide, aes(x =BIN_START/1e6, y = TajimaD), color = "darkorange2") +
  geom_line(data = LAO.tajd.window.wide, aes(x =BIN_START/1e6, y = TajimaD), color = "blue4") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "Tajima's D\n(10 Kbp)") +
  theme_classic(base_size = 12)

p4 <-ggplot() +
  geom_rect(aes(xmin = CK1.start, xmax = CK1.end, ymin = -Inf, ymax = Inf), fill = "forestgreen", alpha = .5) +
  geom_line(data =  STLO.window.het.wide , aes(x =start/1e6, y = mean.hz, color = "GL Odd")) +
  geom_line(data =  LAO.window.het.wide , aes(x =start/1e6, y = mean.hz, color = "BC Odd"))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual( name = "", values = c("BC Odd" = "Blue4", 
                                            "GL Odd" = "darkorange2")) +
  labs(x = paste("Chromosome", chr, "position (Mbp)"), y = "Heterozygosity\n(10 Kbp)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")


all.p <- p1/p2/p3/p4

all.p

ggsave('~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/plot.CK1.SI.pdf', all.p, height  = 8, width = 8.5, units = "in", dpi = 300 )
