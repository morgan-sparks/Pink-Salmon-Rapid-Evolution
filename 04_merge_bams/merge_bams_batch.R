#=============================================================================================================#
# Script created by Morgan Sparks
# This script: uses 2 loops to write out job submission files iterating over samples names and runs for 
# to merge bam files using picard tools to aligned bam files.
# Usage notes: 
#============================================================================================================#

setwd("/scratch/bell/sparks35/GL_Pink_Salmon/scripts/04_merge_bams/")

file_names <- read.table("/scratch/bell/sparks35/GL_Pink_Salmon/data/population_lists/sample_names.txt", sep ="\r", header = FALSE)
file_names <- as.vector(as.character(file_names[,1]))



for(x in c(1:length(file_names))){ # iterate over file names
    x <- noquote(file_names[x])
    
    file_1 <- rbind(
      "#!/bin/bash",
      "",
      paste("#SBATCH --job-name=mergebams_", x, sep = ""),
      "#SBATCH -A standby",
      "#SBATCH -t 4:00:00",
      "#SBATCH -N 1",
      "#SBATCH -n 64",
      "#SBATCH --mail-type=FAIL",
      "#SBATCH --mail-user=sparks35@purdue.edu",
      "",
      "module purge",
      "module load bioinfo",
      "module load picard-tools/2.20.8",
      "module load samtools/1.12",
      "",
      "PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon",
      "ADDRG=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0//03_addedRG_bam/alignment_even",
      "MERGED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/04_merged_bams/alignment_even",
      "",
      "",
      
      paste("java -Xmx124g -jar /group/bioinfo/apps/apps/picard-tools-2.9.0/picard.jar MergeSamFiles \\"),
      paste("I=$ADDRG/", x,"_run1_Ogor1.0_addedRG.bam \\", sep = ""),
      paste("I=$ADDRG/", x,"_run2_Ogor1.0_addedRG.bam \\", sep = ""),
      paste("SORT_ORDER=coordinate \\"),
      paste("O=$MERGED/", x,"_Ogor1.0_merged.bam", sep = "" ),
      "",
      paste("samtools flagstat -@ 126 -O tsv $MERGED/",x,"_Ogor1.0_merged.bam > $MERGED/merged_flagstat/",x, "_rgroups_flagstat_out.txt", sep ="")
     
    )
    
    writeLines(file_1, paste("mergebams_",x,".sh", sep = ""))
  }




