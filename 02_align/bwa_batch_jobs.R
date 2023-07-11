#=============================================================================================================#
# Script created by Morgan Sparks
# This script: uses 2 loops to write out job submission files iterating over samples names and runs for 
# bwa mem
# Usage notes: 
#============================================================================================================#


setwd("/scratch/bell/sparks35/GL_Pink_Salmon/scripts/02_alignment_batch_scripts/alignment_v2")

file_names <- read.table("/scratch/bell/sparks35/GL_Pink_Salmon/data/population_lists/sample_names.txt", sep ="\r", header = FALSE)
file_names <- as.vector(as.character(file_names[,1]))

runs <- c("run1", "run2")

for(z in 1:length(runs)){ #iterate over run1 and run2
  z <- noquote(runs[z])
  for(x in 1:length(file_names)){ # iterate over file names
    x <- noquote(file_names[x])
  file_1 <- rbind(
    "#!/bin/bash",
    "",
    paste("#SBATCH --job-name=align_", x,"_", z, sep = ""),
    "#SBATCH -A standby",
    "#SBATCH -t 4:00:00",
    "#SBATCH -N 1",
    "#SBATCH -n 126",
    "#SBATCH --mail-type=FAIL",
    "#SBATCH --mail-user=sparks35@purdue.edu",
    "",
    "module purge",
    "module load bioinfo",
    "module load bwa/0.7.17",
    "",
    "PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon",
    "",
    "cd $PROJHOME/data/seqs/aligned_reads_Ogor1.0/",
    "",
 paste("bwa mem -t 126 -M $PROJHOME/data/assemblies/Ogor_1.0/GCA_017355495.1_Ogor_1.0_genomic.fna $PROJHOME/data/seqs/trimmed_seqs/trimmed_fastq_v2/trimmed_paired_",x,"_R1_",z,".fastq"," ","$PROJHOME/data/seqs/trimmed_seqs/trimmed_fastq_v2/trimmed_paired_",x,"_R2_",z,".fastq"," ","> $PROJHOME/data/seqs/aligned_reads_Ogor1.0/02_batch_output/alignment_v2/",x,"_",z,"_Ogor1.0_aligned.bam", sep = "")
    )
  
  writeLines(file_1, paste("align_",x,"_",z,".sh", sep = ""))
  }
}



