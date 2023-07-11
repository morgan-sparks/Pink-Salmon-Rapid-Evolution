#=============================================================================================================#
# Script created by Morgan Sparks
# This script loops over chromosome names to apply the GATK hardfilter to the vcf for each chromosome created
# using genotype gvcfs in GATK
#============================================================================================================#


setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon/scripts/09_filterVCFs")

chrom_names <- read.table("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/even_assembly_chromosomes.txt", sep ="\r", header = FALSE)
chrom_names <- as.vector(as.character(chrom_names[,1]))


for(x in c(1:length(chrom_names))){ # iterate over file names
  x <- noquote(chrom_names[x])
  
  file_1 <- rbind(
    "#!/bin/bash",
    "",
    paste("#SBATCH --job-name=hardfilter_", x, sep = ""),
    "#SBATCH -A beagle",
    "#SBATCH -t 24:00:00",
    "#SBATCH -N 1",
    "#SBATCH -n 10",
    "#SBATCH --mail-type=FAIL",
    "#SBATCH --mail-user=sparks35@purdue.edu",
    "",
    "module purge",
    "",
    
    "PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon",
    "CALLSNPS=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/08_genotypeGVCFs",
    "FILTERED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/09_filterVCFs/01_GATK_hardfilter",
    "",
    "",
    "#select SNP variants only",
    "/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx20g \" SelectVariants \\",
    paste("-V $CALLSNPS/allVariants/", x, ".even.vcf.gz \\", sep = ""),
    "-select-type SNP \\",
    paste("-O $CALLSNPS/SNPsOnly/",x, "_snps.vcf.gz ", sep = ""),
    "",

    "#filter using GATK recommended hard filter thresholds https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants",
    "/home/sparks35/gatk-4.2.2.0/gatk --java-options \"-Xmx20g \" VariantFiltration \\",
    paste("-V $CALLSNPS/allSites/",x, ".vcf.gz \\", sep = ""),
    "-filter \"QD < 2.0\" --filter-name \"QD2\" \\",
    "-filter \"SOR > 3.0\" --filter-name \"SOR3\" \\",
    "-filter \"FS > 60.0\" --filter-name \"FS60\" \\",
    "-filter \"MQ < 40.0\" --filter-name \"MQ40\" \\",
    "-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \\",
    "-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \\",
    paste("-O $FILTERED/",x,"_even_GATKfilter_snps.vcf.gz", sep = "")
   
  )
  
  writeLines(file_1, paste("hardfilter_",x,".sh", sep = ""))
}




