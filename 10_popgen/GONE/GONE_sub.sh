#!/bin/bash
#SBATCH --job-name=GONE_itr
#SBATCH -A beagle
#SBATCH -t 300:00:00
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sparks35@purdue.edu

cd /scratch/bell/sparks35/GL_Pink_Salmon/scripts/10_popgen/GONE/GONE_itr/GONE/Linux/

# remove old files so you don't append them
rm -rf *_totalOutput.allloci.txt

for i in {1..500}
do
echo $i
    for FILE in STLO LAO
    do
    bash script_GONE.sh ${FILE}samples_allloci
    tail -n +3 Output_Ne_${FILE}samples_allloci >> ${FILE}_totalOutput.allloci.txt
    done
done


