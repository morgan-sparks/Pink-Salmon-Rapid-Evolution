#!/bin/bash
#SBATCH --job-name=make_genomic_dict
#SBATCH -A beagle
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sparks35@purdue.edu
cd $SLURM_SUBMIT_DIR
module purge


/home/sparks35/gatk-4.2.2.0/gatk CreateSequenceDictionary -R /scratch/bell/sparks35/GL_Pink_Salmon/data/assemblies/OgorEven_v1.0/GCF_021184085.1_OgorEven_v1.0_genomic.fna
