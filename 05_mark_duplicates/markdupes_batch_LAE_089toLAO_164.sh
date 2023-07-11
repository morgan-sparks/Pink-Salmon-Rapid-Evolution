#!/bin/bash

#SBATCH --job-name=markdupes_batch_LAE_089toLAO_170
#SBATCH -A highmem
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 126
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sparks35@purdue.edu

module purge
module load bioinfo

PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
MDUPES=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/05_mark_dupes
MERGED=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/04_merged_bams/alignment_even


#process 1
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAE_089_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAE_089_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAE_089_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=11'

#process 2
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAE_093_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAE_093_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAE_093_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=11'

#process 3
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAE_095_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAE_095_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAE_095_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=11'

#process 4
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAE_098_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAE_098_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAE_098_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=11'

#process 5
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAE_099_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAE_099_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAE_099_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=11'

#process 6
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAE_100_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAE_100_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAE_100_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=11'

#process 7
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAO_161_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAO_161_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAO_161_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=11'

#process 8
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAO_162_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAO_162_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAO_162_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=10'

#process 9
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAO_163_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAO_163_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAO_163_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=10'

#process 10
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAO_164_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAO_164_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAO_164_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=10'

#process 11
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAO_165_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAO_165_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAO_165_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=10'

#process 12
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx80G -Djava.io.tmpdir=/scratch/bell/sparks35/tmpdir" MarkDuplicatesSpark \
-I $MERGED/LAO_170_Ogor1.0_merged.bam \
-O $MDUPES/spark_out/LAO_170_Ogor1.0_dupmarked.bam \
-M $MDUPES/metrics_out/LAO_170_Ogor1.0_dupmarked_metrics.txt \
--conf 'spark.executor.cores=10'

