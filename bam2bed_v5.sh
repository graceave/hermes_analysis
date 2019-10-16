#!/bin/bash                     
#SBATCH --job-name=bam2bed
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-8
#SBATCH --mem=8GB
#SBATCH --mail-type=END
#SBATCH -o b2b.%N.%j.out        # STDOUT
#SBATCH -e b2b.%N.%j.err     # STDERR
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
cd d_$NAME

## convert bam to bed, keeping only read 1 and reads with mapq >5
module purge
module load samtools/intel/1.9
module load bedtools/intel/2.27.1
samtools view -bf 0x40 -q 5 ${NAME}_bktk_sort.bam | bedtools bamtobed -i stdin > ${NAME}_r1.bed

## sort by chromosome and then by start position
sort -k 1,1 -k2,2n ${NAME}_r1.bed > ${NAME}_sorted.bed

