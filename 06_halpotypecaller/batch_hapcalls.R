q#=============================================================================================================#
# Script created by Morgan Sparks
# This script: uses a loop to loop file names to output a bash script that call haplotypecaller
# from GATK. The script loops through chromosomes of pink salmon to call haplotypes in a vcf file
# and then concatenates those into a single file and indexes them at the end of the script with bcftools.
# Usage notes: 
#============================================================================================================#


setwd("/scratch/bell/sparks35/GL_Pink_Salmon/scripts/06_haplotypecaller")

file_names <- read.table("/scratch/bell/sparks35/GL_Pink_Salmon/data/population_lists/sample_names.txt", sep ="\r", header = FALSE)
file_names <- as.vector(as.character(file_names[,1]))



for(x in c(1:length(file_names))){ # iterate over file names
  x <- noquote(file_names[x])
  
  file_1 <- rbind(
    "#!/bin/bash",
    "",
    paste("#SBATCH --job-name=hapcall_", x, sep = ""),
    "#SBATCH -A standby",
    "#SBATCH -t 4:00:00",
    "#SBATCH -N 1",
    "#SBATCH -n 126",
    "#SBATCH --mail-type=FAIL",
    "#SBATCH --mail-user=sparks35@purdue.edu",
    "",
    "module purge",
    "module load bioinfo",
    "module load bcftools/1.11",
    'module load samtools/1.7',
    "",
    
    "PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon",
    "ASSEMBLY=/scratch/bell/sparks35/GL_Pink_Salmon/data/assemblies/OgorEven_v1.0/GCF_021184085.1_OgorEven_v1.0_genomic.fna",
    "MDUPES=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/05_mark_dupes/spark_out",
    "HAPCALLS=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/06_hap_calls",
    "",
    "",
    "while read -a line",
    "do",
    "/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx9g -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir\" HaplotypeCaller \\",
    paste("-I $MDUPES/", x, "_Ogor1.0_dupmarked.bam \\", sep = ""),
    paste("-O $HAPCALLS/${line[0]}_",x, "_OgorEven.g.vcf.gz \\", sep = ""),
    "-R $ASSEMBLY \\",
    "-ERC GVCF \\",
    "-L ${line[0]} &",
    "done < $PROJHOME/data/assemblies/OgorEven_v1.0/even_assembly_chromosomes.txt",
    "wait",
    "",

    "cd /scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/06_hap_calls",
    "",
    "bcftools concat -Oz \\",
    paste("NC_060173.1_",x,"_OgorEven.g.vcf.gz NC_060174.1_",x,"_OgorEven.g.vcf.gz NC_060175.1_",x,"_OgorEven.g.vcf.gz NC_060176.1_",x,"_OgorEven.g.vcf.gz \\", sep =""),
    paste("NC_060177.1_",x,"_OgorEven.g.vcf.gz NC_060178.1_",x,"_OgorEven.g.vcf.gz NC_060179.1_",x,"_OgorEven.g.vcf.gz NC_060180.1_",x,"_OgorEven.g.vcf.gz \\", sep =""),
    paste("NC_060181.1_",x,"_OgorEven.g.vcf.gz NC_060182.1_",x,"_OgorEven.g.vcf.gz NC_060183.1_",x,"_OgorEven.g.vcf.gz NC_060184.1_",x,"_OgorEven.g.vcf.gz \\", sep =""),
    paste("NC_060185.1_",x,"_OgorEven.g.vcf.gz NC_060186.1_",x,"_OgorEven.g.vcf.gz NC_060187.1_",x,"_OgorEven.g.vcf.gz NC_060188.1_",x,"_OgorEven.g.vcf.gz \\", sep =""),
    paste("NC_060189.1_",x,"_OgorEven.g.vcf.gz NC_060190.1_",x,"_OgorEven.g.vcf.gz NC_060191.1_",x,"_OgorEven.g.vcf.gz NC_060192.1_",x,"_OgorEven.g.vcf.gz \\", sep =""),
    paste("NC_060193.1_",x,"_OgorEven.g.vcf.gz NC_060194.1_",x,"_OgorEven.g.vcf.gz NC_060195.1_",x,"_OgorEven.g.vcf.gz NC_060196.1_",x,"_OgorEven.g.vcf.gz \\", sep =""),
    paste("NC_060197.1_",x,"_OgorEven.g.vcf.gz NC_060198.1_",x,"_OgorEven.g.vcf.gz  > ",x,"_OgorEven_hapcalls.g.vcf.gz", sep =""),
    "",
    paste("tabix --csi ",x,"_OgorEven_hapcalls.g.vcf.gz", sep =""),
    "",
    "# remove NC_060199.1_SAMPLE_OgorEven.g.vcf.gz because it is Y chrom and will fail if not present"
  )
  
  writeLines(file_1, paste("hapcall_",x,".sh", sep = ""))
}



