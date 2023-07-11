#!/bin/bash
#SBATCH --job-name=eigenGWAS_LAOvSTLO
#SBATCH -A beagle
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -t 330:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sparks35@purdue.edu

GWASDIR=/scratch/bell/sparks35/GL_Pink_Salmon/data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/eigenGWAS

java -jar ~/GEAR/gear.jar eigengwas --bfile $GWASDIR/LAO_STLO_samples --ev 2 --thread-num 15 --out $GWASDIR/LAO_STLO_samples

