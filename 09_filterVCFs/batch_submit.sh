while read -a line
do 
sbatch ${line[0]}
done < job_names.txt