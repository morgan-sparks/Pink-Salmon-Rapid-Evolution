library(rtracklayer); library(dplyr)
setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon")

dat <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/LAOvSTLOoutlier_windows_all_individuals.txt", header = TRUE)
#dat <- read.table("./scripts/10_popgen/Fst/outlier_windows_GL2vs3yrs.txt", header = TRUE)

dat$c <- dplyr::recode(as.character(dat$c),
                       "1"  = "NC_060173.1" ,
                       "2"  = "NC_060174.1" ,
                       "3"  = "NC_060175.1" ,
                       "4"  = "NC_060176.1" ,
                       "5"  = "NC_060177.1" ,
                       "6"  = "NC_060178.1" ,
                       "7"  = "NC_060179.1" ,
                       "8"  = "NC_060180.1" ,
                       "9"  = "NC_060181.1" ,
                       "10" ="NC_060182.1",
                       "11" ="NC_060183.1",
                       "12" ="NC_060184.1",
                       "13" ="NC_060185.1",
                       "14" ="NC_060186.1",
                       "15" ="NC_060187.1",
                       "16" ="NC_060188.1",
                       "17" ="NC_060189.1",
                       "18" ="NC_060190.1",
                       "19" ="NC_060191.1",
                       "20" ="NC_060192.1",
                       "21" ="NC_060193.1",
                       "22" ="NC_060194.1",
                       "23" ="NC_060195.1",
                       "24" ="NC_060196.1",
                       "25" ="NC_060197.1",
                       "26" ="NC_060198.1")



# load gff
gff <- ape::read.gff("./data/assemblies/OgorEven_v1.0/GCF_021184085.1_OgorEven_v1.0_genomic.gff")
# remove regions (e.g. chromomosomes and scaffolds)
gff1 <- gff[-which(gff$type == "region"),]


OUT5 = NULL
for (i in 1:nrow(dat)){
  outlier <- dat[i,]
  
  #subset gff to match only chromosome of interest
  gff_small <- gff1[which(gff1$seqid == outlier$c),]
  
  #gff if window is fully between start and stop of an attribute
  full_match <- gff_small[which(gff_small$start >= outlier$real.starts & gff_small$end <= outlier$ends),]

  #full <- cbind(cbind(rep(outlier, 2)),  match = rep("full", nrow(full_match)), full_match)

  
  #gff if window is starts before start but ends before end
 
  start_overlap <- gff_small[which(gff_small$start <= outlier$real.starts & gff_small$ends <= outlier$real.starts& gff_small$end <= outlier$ends),]
  
  #gff if window is starts after start but runs over end
  end_overlap <- gff_small[which(gff_small$start >= outlier$real.starts & gff_small$start <= outlier$ends & gff_small$end >= outlier$ends),]
  
  all_results <- rbind(full_match, start_overlap, end_overlap)
  
  outlier_all <-cbind(outlier_num = rep(i, nrow(all_results)), 
              outlier_start = rep(outlier[,2], nrow(all_results)),
              outlier_end = rep(outlier[,3], nrow(all_results)),
              zfst = rep(outlier[,4], nrow(all_results)),
              match = c(rep("full", nrow(full_match)), rep("start", nrow(start_overlap)), rep("end", nrow(end_overlap))),
              all_results)
  
  OUT5 <-rbind(OUT5, outlier_all)
  
  
}

outlier_hits <- OUT5

### write out outlier hits

write.table(outlier_hits, "./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers/LAOvSTLO_outlier_gffhits.csv", sep = ",", col.names=TRUE, row.names=FALSE, append=FALSE)

####
outlier_hits[order(-outlier_hits$zfst),]

dat[order(-dat$z.fst),]

genes<- outlier_hits[outlier_hits$type == "gene",]


genes %>% 
  arrange(desc(zfst)) %>% 
  select(attributes) %>% 
  top_n(10)

outlier_hits %>% 
  filter(outlier_num == 43)


##### chr10 inv genes

chr10_genes <- gff[gff$seqid == "NC_060182.1" & gff$type == "gene"  & gff$start >= 10e6 & gff$end <= 35e6,]


