library(ggplot2); library(patchwork)
setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon/scripts/10_popgen/misc")

list.files()

file_names <- read.table("./combinedBamsFileName.txt", sep = "\n", header = FALSE, stringsAsFactors = FALSE)

samples <- unlist(strsplit(file_names$V1, "_Ogor"))

samples <- samples[seq(from = 1, to = length(samples), by = 2)]

Q20_reads <- read.table("./combinedBamsQ20Reads.txt", sep = "\n", header = FALSE, stringsAsFactors = FALSE)

Q60_reads <- read.table("./combinedBamsQ60Reads.txt", sep = "\n", header = FALSE, stringsAsFactors = FALSE)

mapped_reads <- read.table("./combinedBamsMappedReads.txt", sep = "\n", header = FALSE, stringsAsFactors = FALSE)

total_reads <- read.table("./combinedBamsTotalReads.txt", sep = "\n", header = FALSE, stringsAsFactors = FALSE)


read_data <- cbind(samples, Q20_reads$V1, Q60_reads$V1, mapped_reads$V1, total_reads$V1)

colnames(read_data) <- c("samples", "Q20", "Q60", "mapped_reads", "total_reads")


read_data <- data.frame(read_data)

for( i in c(2:5)){
read_data[,i] <- as.numeric(as.character(read_data[,i]))
}

Q20_plot <- ggplot(data = read_data)+ 
  geom_bar(aes(x= samples, y = ((Q20/total_reads)*100)), stat = "identity") +
  theme_classic() +
  labs(x = "Sample", y = "Percent >/= MapQ 20", title = paste("Mean =", round(mean((read_data$Q20/read_data$total_reads)*100), 2), "%", sep = " ")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Q60_plot <- ggplot(data = read_data)+ 
  geom_bar(aes(x= samples, y = ((Q60/total_reads)*100)), stat = "identity") +
  theme_classic() +
  labs(x = "Sample", y = "Percent >/= MapQ 60", title = paste("Mean =", round(mean((read_data$Q60/read_data$total_reads)*100), 2), "%", sep = " ")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mapread_plot <- ggplot(data = read_data)+ 
  geom_bar(aes(x= samples, y = (mapped_reads/total_reads)*100), stat = "identity") +
  theme_classic() +
  labs(x = "Sample", y = "Percent Mapped Reads", title = paste("Mean =", round(mean((read_data$mapped_reads/read_data$total_reads)*100), 2), "%", sep = " ")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

totalread_plot<- ggplot(data = read_data)+ 
  geom_bar(aes(x= samples, y = total_reads), stat = "identity") +
  theme_classic() +
  labs(x = "Sample", y = "Mapped reads", title ="Total mapped reads") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

(Q20_plot + Q60_plot) / (mapread_plot + totalread_plot)


ggplot(data = read_data)+ 
  geom_bar(aes(x= samples, y = (total_reads - mapped_reads)), stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#######
# plot per sample avg mapq w/ stdev

avg_MapQ <- read.table("./BamAvgMapQ.txt", sep = "\n", header = FALSE, stringsAsFactors = FALSE)

avg_MapQ < cbind(samples = samples, avg_MapQ = avg_MapQ)  
  
  