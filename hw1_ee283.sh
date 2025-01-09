#!/bin/bash
#SBATCH --job-name=test    ## Name of the job.
#SBATCH -A CLASS-ECOEVO283       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --cpus-per-task=1  ## number of cores the job needs
#SBATCH --error=hw1_ee283.%J.err
#SBATCH --output=hw1_ee283.%J.out

wget https://wfitch.bio.uci.edu/~tdlong/problem1.tar.gz
tar -xvf problem1.tar.gz
rm problem1.tar.gz

head problem1/p.txt>1p.txt
tail -1 1p.txt 
rm 1p.txt

head problem1/f.txt>1f.txt
tail -1 1f.txt
rm 1f.txt

sleep 2m	# wait 2 minutes
