setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon")

library(GenotypePlot); library(tidyverse)
vcf <- vcfR::read.vcfR("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/GL_Pink_Salmon_chr10.recode.vcf")

# vcf.mat <- data.frame(vcf@fix)
# vcf.mat$POS <- as.numeric(vcf.mat$POS)
# vcf.mat <- vcf.mat[vcf.mat$POS <= 40e6,]
# tail(vcf.mat) # the 87332 row is the last before 40 Mbp



set.seed(8675309)
rand.pos <- sort(sample(1:87332, size = 5000, replace = FALSE))

my_vcf <- vcf[rand.pos,]

rm(vcf)


inds <- data.frame(colnames(my_vcf@gt)[2:length(colnames(my_vcf@gt))])
colnames(inds) <- "samples"

pop <- inds %>% 
  separate(samples, sep = "_", into = c("pop", "ind")) %>% 
  dplyr::select("pop")

our_popmap <- data.frame(cbind(ind = inds, pop = pop))
colnames(our_popmap)[1] <- "ind"

our_popmap$pop <- recode(our_popmap$pop, "LAO" = "BC Odd",
                                     "LAE" = "BC Even",
                                     "STLO" = "GL Odd",
                                     "STLO3" = "GL Odd 3",
                                     "STLE" = "GL Even"
                         )
### pop maps for pixy

# LAOvSTLO.popmap <- our_popmap[our_popmap$pop %in% c("BC Odd", "GL Odd"),]
# STLOvSTLE.popmap  <- our_popmap[our_popmap$pop %in% c("GL Odd", "GL Even"),]
# # 
# write.table("./scripts/10_popgen/pixy/LAOvSTLO.popmap.txt", x = LAOvSTLO.popmap, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
# 
# write.table("./scripts/10_popgen/pixy/STLOvSTLE.popmap.txt", x = STLOvSTLE.popmap, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

X <- cbind(pops = our_popmap, X)

X.long <- pivot_longer(X, cols= c(3:ncol(X)), names_to = "pos", values_to = "gntyp")
X.long$pos <- gsub(pattern = "X",  replacement = "", X.long$pos)
X.long$pos <- as.numeric(X.long$pos)
X.long$gntyp<- as.character(X.long$gntyp)
X.long <- na.omit(X.long)



gntyp.pos <- unique(sort(X.long$pos))



gntyp.breaks <- c(gntyp.pos[1], 
                  gntyp.pos[length(gntyp.pos)])

gntyp.breaks <- as.character(gntyp.breaks)



gntyp.plot <- ggplot() +
  geom_tile(data = X.long, aes(x = as.character(pos), y=pops.ind, fill= gntyp)) +
  scale_fill_manual(values = gntyp.pal,
                    labels = c("HOM Ref", "HET", "HOM Alt")) +
  scale_x_discrete(breaks = gntyp.breaks, 
                   labels = c("0", "55")) +
  geom_hline(yintercept = c(29, 59, 89, 119)) +
  labs(x = "Position (Mbp)", y = "") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank())

ggsave("~/Downloads/gntyp.plot.png", gntyp.plot, height = 4, width = 12, units = "in", dpi = 300)

old.pal <- c("#C92D59","#FCD225","#300060")
gntyp.pal<-c("#faf3dd","#a5be00", "#300060")


new_plot <- genotype_plot(vcf_object  =  my_vcf,                                      # chr or scaffold ID
                          popmap = our_popmap,                              # population membership
                          cluster        = FALSE,                           # whether to organise haplotypes by PCA clustering
                          snp_label_size = 10e6,                          # breaks for position labels, eg. plot a position every 100,000 bp
                          colour_scheme=gntyp.pal,                       # character vector of colour values
                          invariant_filter = TRUE)                      # Filter any invariant sites before plotting
              

geno.plot <- new_plot$genotypes +
  #annotate("rect", xmin =6403689, xmax = 35403689, ymin = -Inf, ymax = Inf, fill = "blue4", alpha = 0.25 ) 
  #geom_rect(aes(xmin= 6403689, xmax = 35403689, ymin = -Inf, ymax = Inf), color = "white")
  geom_vline(xintercept = c(6403689,35403689), color = "black")

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/CH 3/Figures/Fig. 2b.png",geno.plot, height = 3.5, width = 17, units = "in", dpi =1200)
ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/CH 3/Figures/Fig. 2b label.png", new_plot$positions, height = 2, width = 8.5, units = "in", dpi =1200)





png("~/Downloads/chr10_gntyp2s.png",  height = 10, width =30, unit = "in",res = 300)
new_plot$genotypes +
  geom_vline(xintercept = 6403689)
dev.off()
 


# snps <-nrow(my_vcf)
# new_plot2 <- genotype_plot(vcf_object  =  my_vcf[1:snps/2,],                                      # chr or scaffold ID
#                           
#                           popmap = our_popmap,                              # population membership
#                           cluster        = FALSE,                           # whether to organise haplotypes by PCA clustering
#                           snp_label_size = 5e6,                          # breaks for position labels, eg. plot a position every 100,000 bp
#                           colour_scheme=c("#faf3dd","#a5be00", "#513b56"),   # character vector of colour values
#                           invariant_filter = TRUE) 

### try haplotype 

vcf <- vcfR::read.vcfR("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/inversion/EHH/STLOsamples_bypop_chr10_phased.vcf")

my_vcf <- vcf[1:69961,]

inds$pop <- "GL Odd"

our_popmap <- inds


p1 <- genotype_plot(vcf_object  =  my_vcf,                                       # chr or scaffold ID
              plot_phased = TRUE,
              popmap = our_popmap,                              # population membership
              cluster        = FALSE,                           # whether to organise haplotypes by PCA clustering
              snp_label_size = 10e6,                          # breaks for position labels, eg. plot a position every 100,000 bp
              #colour_scheme=c("dodgerblue","darkred","blue2"),   # character vector of colour values
              colour_scheme=c("#FCD225","#300060"),
              invariant_filter = TRUE)  


vcf.2 <- vcfR::read.vcfR("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/inversion/EHH/LAOsamples_bypop_chr10_phased.vcf")

my_vcf.2 <- vcf.2[1:69961,]

inds.2 <- data.frame(colnames(my_vcf.2@gt)[2:length(colnames(my_vcf.2@gt))])
inds.2$pop <- "BC Odd"

our_popmap.2 <- inds.2


p2 <- genotype_plot(vcf_object  =  my_vcf.2,                                       # chr or scaffold ID
              plot_phased = TRUE,
              popmap = our_popmap.2,                              # population membership
              cluster        = FALSE,                           # whether to organise haplotypes by PCA clustering
              snp_label_size = 10e6,                          # breaks for position labels, eg. plot a position every 100,000 bp
              #colour_scheme=c("dodgerblue","darkred","blue2"),   # character vector of colour values
              colour_scheme=c("#FCD225","#300060"),
              invariant_filter = TRUE)

library(patchwork)

haplo.genos <- p2$genotypes / p1$genotypes

ggsave("~/Downloads/haplo.genotypes.png", haplo.genos, width  = 11, height = 6, units = "in", dpi = 600)
