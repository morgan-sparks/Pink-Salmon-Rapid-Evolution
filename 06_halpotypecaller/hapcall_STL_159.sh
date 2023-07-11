#!/bin/bash

#SBATCH --job-name=hapcall_STL_159
#SBATCH -A standby
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 126
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge
module load bioinfo
module load bcftools/1.11
module load samtools/1.7

PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
ASSEMBLY=/scratch/bell/sparks35/GL_Pink_Salmon/data/assemblies/OgorEven_v1.0/GCF_021184085.1_OgorEven_v1.0_genomic.fna
MDUPES=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/05_mark_dupes/spark_out
HAPCALLS=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/06_hap_calls


while read -a line
do
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx9g -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" HaplotypeCaller \
-I $MDUPES/STL_159_Ogor1.0_dupmarked.bam \
-O $HAPCALLS/${line[0]}_STL_159_OgorEven.g.vcf.gz \
-R $ASSEMBLY \
-ERC GVCF \
-L ${line[0]} &
done < $PROJHOME/data/assemblies/OgorEven_v1.0/even_assembly_chromosomes.txt
wait

cd /scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/06_hap_calls

bcftools concat -Oz \
NC_060173.1_STL_159_OgorEven.g.vcf.gz NC_060174.1_STL_159_OgorEven.g.vcf.gz NC_060175.1_STL_159_OgorEven.g.vcf.gz NC_060176.1_STL_159_OgorEven.g.vcf.gz \
NC_060177.1_STL_159_OgorEven.g.vcf.gz NC_060178.1_STL_159_OgorEven.g.vcf.gz NC_060179.1_STL_159_OgorEven.g.vcf.gz NC_060180.1_STL_159_OgorEven.g.vcf.gz \
NC_060181.1_STL_159_OgorEven.g.vcf.gz NC_060182.1_STL_159_OgorEven.g.vcf.gz NC_060183.1_STL_159_OgorEven.g.vcf.gz NC_060184.1_STL_159_OgorEven.g.vcf.gz \
NC_060185.1_STL_159_OgorEven.g.vcf.gz NC_060186.1_STL_159_OgorEven.g.vcf.gz NC_060187.1_STL_159_OgorEven.g.vcf.gz NC_060188.1_STL_159_OgorEven.g.vcf.gz \
NC_060189.1_STL_159_OgorEven.g.vcf.gz NC_060190.1_STL_159_OgorEven.g.vcf.gz NC_060191.1_STL_159_OgorEven.g.vcf.gz NC_060192.1_STL_159_OgorEven.g.vcf.gz \
NC_060193.1_STL_159_OgorEven.g.vcf.gz NC_060194.1_STL_159_OgorEven.g.vcf.gz NC_060195.1_STL_159_OgorEven.g.vcf.gz NC_060196.1_STL_159_OgorEven.g.vcf.gz \
NC_060197.1_STL_159_OgorEven.g.vcf.gz NC_060198.1_STL_159_OgorEven.g.vcf.gz  > STL_159_OgorEven_hapcalls.g.vcf.gz

tabix --csi STL_159_OgorEven_hapcalls.g.vcf.gz

# remove NC_060199.1_SAMPLE_OgorEven.g.vcf.gz because it is Y chrom and will fail if not present
