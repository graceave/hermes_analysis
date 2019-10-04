#!/bin/bash                     
#SBATCH --job-name=5dedup
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0
#SBATCH --mem=8GB
#SBATCH --mail-type=END
#SBATCH -o tn_dedup.%N.%j.out        # STDOUT
#SBATCH -e tn_dedup.%N.%j.err     # STDERR
#SBATCH --mail-user=ga824@nyu.edu
##SBATCH --dependency=afterok:

##set variables
REF="/genomics/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GAP1/GCF_000146045.2_R64_genomic_GAP1.fna"

##capture the output of a command line and store it in a variable
ALL_R1=($(ls /scratch/cgsb/gencore/out/Gresham/2019-09-17_HTFWFAFXY/merged/*n01*))

##this is going to assign the variables to file names
FORWARD=${ALL_R1[$SLURM_ARRAY_TASK_ID]}
NAME=${FORWARD:76:-22}

echo $NAME

## change into directory for this specific sample
cd d_${NAME}


## filter out multimapping reads (mq<5) and deduplicate using umis
module purge
module load umi_tools/0.5.5
umi_tools dedup --extract-umi-method=read_id --mapping-quality=5 --paired --output-stats=${NAME} --method=directional --stdin=${NAME}_bktk_sort.bam --log=${NAME}_dedup.log > ${NAME}_dedup.bam 
