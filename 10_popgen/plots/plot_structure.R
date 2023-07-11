# Header ------------------------------------------------------------------
# structure plot using fastSTRUCUTRE, ran K=1-5 for all GL Pink salmon samples.
# It determined K = 3 was the best option.

library(tidyverse)

setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon")

### data

struct <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/structure/results/test.3.meanQ")
pops <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/structure/data/GL_PinkSalmon_jointfiltered_10KrandomSNPS.fam")

### set up data

allpops.struct <- data.frame(cbind(pops[,1:2], struct))
colnames(allpops.struct) <- c("pop", "sample", "GL", "BC Even", "BC Odd")

allpops.struct <- allpops.struct %>%  unite(allpops.struct, pop, sample, sep = "_") 
colnames(allpops.struct)[1] <- "sample"

allpops.struct <- allpops.struct %>% pivot_longer(cols =2:4, names_to = "assignment", values_to = "value")


p.struct <- ggplot(allpops.struct, aes(x = sample, y =value*100, fill = assignment)) +
  geom_bar(position="stack", stat="identity",width = 1) +
  scale_fill_manual(name = "Population", 
                    values = c("GL" = "darkorange2",
                               "BC Odd" = "blue4",
                               "BC Even" = "darkgreen")) +
  labs(x = "Sample", y = "Assignment probability") +
  coord_cartesian(x = c(0,100)) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x=element_text(angle = -90, hjust=1, size = 10))

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/structure.pdf", p.struct, height = 4, width = 11, units = "in", dpi = 600)
