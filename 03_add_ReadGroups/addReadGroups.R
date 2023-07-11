#=============================================================================================================#
# Script created by Morgan Sparks
# This script: uses 2 loops to write out job submission files iterating over samples names and runs for 
# to add Read Groups using picard tools to aligned bam files.
# Usage notes: 
#============================================================================================================#


setwd("/scratch/bell/sparks35/GL_Pink_Salmon/scripts/03_add_ReadGroups")

file_names <- read.table("/scratch/bell/sparks35/GL_Pink_Salmon/data/population_lists/sample_names.txt", sep ="\r", header = FALSE)
file_names <- as.vector(as.character(file_names[,1]))

runs <- c("run1", "run2")

for(z in 1:length(runs)){ #iterate over run1 and run2
  z <- noquote(runs[z])
 
   for(x in c(1:length(file_names))){ # iterate over file names
    x <- noquote(file_names[x])
    
    file_1 <- rbind(
      "#!/bin/bash",
      "",
      paste("#SBATCH --job-name=addreadgroups_", x,"_", z, sep = ""),
      "#SBATCH -A standby",
      "#SBATCH -t 4:00:00",
      "#SBATCH -N 1",
      "#SBATCH -n 20",
      "#SBATCH --mail-type=FAIL",
      "#SBATCH --mail-user=sparks35@purdue.edu",
      "",
      "module purge",
      "module load bioinfo",
      "module load picard-tools/2.20.8",
      "",
      "PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon",
      "ALIGNED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/02_aligned_reads/alignment_even",
      "RGADD=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/03_addedRG_bam/alignment_even",
      "",
      "",
      
      paste("java -Xmx40g -jar /group/bioinfo/apps/apps/picard-tools-2.9.0/picard.jar AddOrReplaceReadGroups \\"),
      paste("I=$ALIGNED/", x,"_", z,"_Ogor1.0_aligned.bam \\", sep = ""),
      paste("O=$RGADD/", x,"_", z,"_Ogor1.0_addedRG.bam \\", sep = "" ),
      paste("SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=", x," RGPU=", z, sep = "")
       )
    
    writeLines(file_1, paste("addRG_",x,"_",z,".sh", sep = ""))
  }
}



