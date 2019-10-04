#!/bin/bash                     
#SBATCH --job-name=5fqc
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --array=0-8
#SBATCH --mem=8GB
#SBATCH --mail-type=END
#SBATCH -o fqc.%N.%j.out        # STDOUT
#SBATCH -e fqc.%N.%j.err     # STDERR
#SBATCH --mail-user=ga824@nyu.edu
##SBATCH --dependency=afterok:

##capture the output of a command line and store it in a variable
ALL_R1=($(ls /scratch/cgsb/gencore/out/Gresham/2019-09-17_HTFWFAFXY/merged/*n01*))

##this is going to assign the variables to file names
FORWARD=${ALL_R1[$SLURM_ARRAY_TASK_ID]}
NAME=${FORWARD:76:-22}

echo ${NAME}

## change to directory for this specific sample
cd d_${NAME}

 ## fastqc of cleaned reads (output of clean_reads_v2.sh)
module purge
module load fastqc/0.11.8
fastqc -t 3 ${NAME}_clean_r1.fastq.gz

