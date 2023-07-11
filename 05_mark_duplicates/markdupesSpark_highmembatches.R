#=============================================================================================================#
# Script created by Morgan Sparks
# This script: uses 2 loops to write out job submission files iterating over samples names and runs for 
# to add Read Groups using picard tools to aligned bam files.
# Usage notes: 
#============================================================================================================#


setwd("/scratch/bell/sparks35/GL_Pink_Salmon/scripts/05_mark_duplicates")

file_names <- read.table("/scratch/bell/sparks35/GL_Pink_Salmon/data/population_lists/sample_names.txt", sep ="\r", header = FALSE)
file_names <- as.vector(as.character(file_names[,1]))

x<-1

# there are 134 files, so this makes batch files in batches of 12 until 132
# and then a final batch or just two

while(x < 134){ # iterate over file names
  
  if(x<132){
  a <- file_names[x]
  b <- file_names[(x+1)]
  c <- file_names[(x+2)]
  d <- file_names[(x+3)]
  e <- file_names[(x+4)]
  f <- file_names[(x+5)]
  g <- file_names[(x+6)]
  h <- file_names[(x+7)]
  i <- file_names[(x+8)]
  j <- file_names[(x+9)]
  k <- file_names[(x+10)]
  l <- file_names[(x+11)]
  
  
  file_1 <- rbind(
    "#!/bin/bash",
    "",
    paste("#SBATCH --job-name=markdupes_batch_", a, "to",l,sep = ""),
    "#SBATCH -A highmem",
    "#SBATCH -t 24:00:00",
    "#SBATCH -N 1",
    "#SBATCH -n 126",
    "#SBATCH --mail-type=FAIL,END",
    "#SBATCH --mail-user=sparks35@purdue.edu",
    "",
    "module purge",
    "module load bioinfo",
    
    "",
    "PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon",
    "MDUPES=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/05_mark_dupes",
    "MERGED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/04_merged_bams/alignment_even",
    "",
    "",
    "#process 1",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",a,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",a,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",a,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=11'",
    "",
    
    "#process 2",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",b,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",b,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",b,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=11'",
    "",
   
    "#process 3",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",c,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",c,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",c,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=11'",
    "",
    "#process 4",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",d,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",d,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",d,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=11'",
    "",
    
    "#process 5",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",e,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",e,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",e,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=11'",
    "",
    
    "#process 6",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",f,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",f,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",f,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=11'",
    "",
    
    "#process 7",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",g,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",g,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",g,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=11'",
    "",
    
    "#process 8",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",h,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",h,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",h,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=10'",
    "",
    
    "#process 9",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",i,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",i,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",i,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=10'",
    "",
    
    "#process 10",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",j,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",j,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",j,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=10'",
    "",
    
    "#process 11",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",k,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",k,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",k,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=10'",
    "",
    
    "#process 12",
    paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
    paste("-I $MERGED/",l,"_Ogor1.0_merged.bam \\", sep = ""),
    paste("-O $MDUPES/spark_out/",l,"_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-M $MDUPES/metrics_out/",l,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
    "--conf 'spark.executor.cores=10'",
    ""
    
   
  )
  
  writeLines(file_1, paste("markdupes_batch_",a,"to",j,".sh", sep = ""))
  x <- x + 12
  }else{
    a <- file_names[x]
    b <- file_names[(x+1)]
    
    
    file_1 <- rbind(
      "#!/bin/bash",
      "",
      paste("#SBATCH --job-name=markdupes_batch_", a, "to",b,sep = ""),
      "#SBATCH -A highmem",
      "#SBATCH -t 24:00:00",
      "#SBATCH -N 1",
      "#SBATCH -n 63",
      "#SBATCH --mail-type=FAIL,END",
      "#SBATCH --mail-user=sparks35@purdue.edu",
      "",
      "module purge",
      "module load bioinfo",
      
      "",
      "PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon",
      "MDUPES=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/05_mark_dupes",
      "MERGED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/04_merged_bams/alignment_even",
      "",
      "",
      
      "#process 1",
      paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx225G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
      paste("-I $MERGED/",a,"_Ogor1.0_merged.bam \\", sep = ""),
      paste("-O $MDUPES/spark_out/",a,"_Ogor1.0_dupmarked.bam \\", sep = ""),
      paste("-M $MDUPES/metrics_out/",a,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
      "--conf 'spark.executor.cores=30'",
      "",
      
      "#process 2",
      paste("/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx225G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" MarkDuplicatesSpark \\"),
      paste("-I $MERGED/",b,"_Ogor1.0_merged.bam \\", sep = ""),
      paste("-O $MDUPES/spark_out/",b,"_Ogor1.0_dupmarked.bam \\", sep = ""),
      paste("-M $MDUPES/metrics_out/",b,"_Ogor1.0_dupmarked_metrics.txt \\", sep = ""),
      "--conf 'spark.executor.cores=30'",
      ""
     
      
    )
    
    writeLines(file_1, paste("markdupes_batch_",a,"to",b,".sh", sep = ""))
    x <- x + 2
    
  }
}





