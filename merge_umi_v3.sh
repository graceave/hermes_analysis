#!/bin/bash                     
#SBATCH --job-name=2mergeumi
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1
#SBATCH --mem=40GB
#SBATCH --mail-type=END
#SBATCH -o umerge2.%N.%j.out        # STDOUT
#SBATCH -e umerge2.%N.%j.err     # STDERR
#SBATCH --mail-user=ga824@nyu.edu
##SBATCH --dependency=afterok:

##set variables
REF="/genomics/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GAP1/GCF_000146045.2_R64_genomic_GAP1.fna"

##capture the output of a command line and store it in a variable
ALL_R1=($(ls /scratch/cgsb/gencore/out/Gresham/2019-09-17_HTFWFAFXY/merged/*n01*))

##this is going to assign the variables to file names
##FORWARD=${ALL_R1[$SLURM_ARRAY_TASK_ID]}
##NAME=${FORWARD:76:-22}
NAME='1657_1'

echo $NAME

## change into directory for this specific sample
cd d_${NAME}

## merge umis with my r script
module purge 
module load r/intel/3.4.2
Rscript /scratch/ga824/hermes_analysis/merge_umi_v3.R ${NAME}_sorted.bed ${SLURM_ARRAY_TASK_ID} ${NAME}_chr${SLURM_ARRAY_TASK_ID}_merge.bed


