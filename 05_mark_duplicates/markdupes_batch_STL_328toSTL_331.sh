#!/bin/bash

#SBATCH --job-name=markdupes_batch_STL_328toSTL_331
#SBATCH -A highmem
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 63
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sparks35@purdue.edu

module purge
module load bioinfo

PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
MDUPES=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/05_mark_dupes
MERGED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/04_merged_bams/alignment_even


#process 1
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx225G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/STL_328_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/STL_328_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/STL_328_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=30'

#process 2
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx225G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/STL_331_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/STL_331_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/STL_331_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=30'

