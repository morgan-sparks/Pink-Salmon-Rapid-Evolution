#!/bin/bash
#SBATCH --job-name=GDBImp_allSamps_ch22_27
#SBATCH -A highmem
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sparks35@purdue.edu

module purge
PROJHOME=/scratch/bell/sparks35/GL_Pink_Salmon
HAPCALLS=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_Ogor1.0/06_hap_calls

while read -a line
do
/home/sparks35/gatk-4.2.2.0/gatk --java-options "-Xmx128g -Xms128g" GenomicsDBImport \
-V $HAPCALLS/${line[0]}_LAO_161_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_162_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_163_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_164_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_165_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_170_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_171_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_173_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_177_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_182_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_183_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_185_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_186_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_187_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_188_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_189_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_191_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_192_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_193_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_194_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_195_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_198_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_199_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_201_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_203_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_204_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_205_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_206_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAO_207_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_101_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_105_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_112_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_113_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_125_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_136_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_159_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_160_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_187_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_189_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_191_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_201_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_205_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_206_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_210_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_226_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_228_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_244_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_247_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_259_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_264_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_270_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_271_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_290_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_304_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_315_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_316_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_318_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_328_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_331_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_006_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_012_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_024_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_030_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_036_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_042_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_053_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_056_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_057_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_058_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_059_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_062_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_064_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_065_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_071_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_076_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_077_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_080_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_081_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_082_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_083_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_086_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_087_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_088_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_089_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_093_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_095_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_098_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_099_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_LAE_100_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_003_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_006_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_007_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_008_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_009_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_011_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_012_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_013_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_014_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_015_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_016_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_017_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_018_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_019_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_020_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_021_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_022_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_024_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_025_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_026_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_028_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_030_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_031_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_033_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_035_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_036_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_038_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_039_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_046_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_048_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_104_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_110_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_111_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_134_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_139_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_153_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_175_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_179_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_194_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_203_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_252_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_302_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_308_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_310_Ogor1.0.g.vcf.gz \
-V $HAPCALLS/${line[0]}_STL_319_Ogor1.0.g.vcf.gz \
--genomicsdb-workspace-path $PROJHOME/data/seqs/aligned_reads_Ogor1.0/07_genomicsDB/all_samps/${line[0]}_genomicsDatabase \
-L ${line[0]} \
--batch-size 30 \
--tmp-dir tmpdir \
--genomicsdb-shared-posixfs-optimizations true \
--reader-threads 4 &
  done < $PROJHOME/data/seqs/aligned_reads_Ogor1.0/chromosome_list_22_27.txt
wait