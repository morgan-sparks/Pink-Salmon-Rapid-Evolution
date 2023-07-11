# How to use GONE program

## Download GONE

First download from GitHub (I did in my home directory on Bell)

`git clone https://github.com/esrud/GONE`

To use in Linux, change permission on executables

`chmod u+x ./GONE/Linux/PROGRAMMES`

## How to run GONE

Gone is easiest to run by just copying the original directory and rerunning with a new copy every time you want to do so. See discussion in https://www.r-bloggers.com/2021/12/estimating-recent-population-history-from-linkage-disequilibrium-with-gone-and-snep/

I copied GONE into my `/scripts` folder in my project directory and made copies for the two populations I wanted to run on `GONE_STLO` AND `GONE_LAO`

The program works by calling .ped and .map files that are stored in the `/GONE/Linux/` directory

First I had to convert my vcf file to plink files using the command. I tried doing this with vcftools but the outputted .ped file had old genotype markers that the program errors out with.

`plink --vcf file.vcf --recode --out file`

Because I am using randomly subsampled vcf files the chromosomes are all out of order and need to be sorted in the .bim file. Resort into a corrected file with:

`plink --bfile file --make-bed --out file.corr`

Finally, plink converts the vcf into a .bed file and not a .ped file, which is necessary for the GONE program so you need to convert that with:

`plink --bfile file.corr --recode tab --out file.corr`

Now the files are ready. To run gone on the cluster you need to execute the `script_GONE.sh` shell script. To do so I created a job script called `GONE_sub.sh` It looks like:

```
#!/bin/bash
#SBATCH --job-name=GONE
#SBATCH -A beagle
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sparks35@purdue.edu
cd $SLURM_SUBMIT_DIR
module purge

bash script_GONE.sh file.corr

```

The key to running this is the `script_GONE.sh` needs to call the file prefix for your corrected .bed and .map  (e.g., file.corr in this case). It should then run quickly and output its files into the `/GONE/Linux/` directory.