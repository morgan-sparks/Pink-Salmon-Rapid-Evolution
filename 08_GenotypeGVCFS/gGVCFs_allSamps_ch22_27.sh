#!/bin/bash
#SBATCH --job-name=gGVCFs_allSamps_ch22_27
#SBATCH -A mcclintock
#SBATCH -N 1
#SBATCH -n 80
#SBATCH -t 330:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge
PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
ASSEMBLY=/scratch/bell/sparks35/GL_Pink_Salmon/data/assemblies/OgorEven_v1.0/GCF_021184085.1_OgorEven_v1.0_genomic.fna
GENOMICSDB=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/07_genomicsDB/all_samps
OUT=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/08_genotypeGVCFs

while read -a line
do
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx10g -Xms10g" GenotypeGVCFs \
-R $ASSEMBLY \
-V gendb://$GENOMICSDB/${line[0]}_genomicsDatabase \
-L ${line[0]} \
-O $OUT/${line[0]}.even.vcf.gz &
done < $PROJHOME/data/seqs/aligned_reads_OgorEven_v1.0/chromosome_list_22_27.txt
wait
