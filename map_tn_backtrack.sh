#!/bin/bash                     
#SBATCH --job-name=v4bt_map
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --array=0
#SBATCH --mem=8GB
#SBATCH --mail-type=END
#SBATCH -o bt_mapv4.%N.%j.out        # STDOUT
#SBATCH -e bt_mapv4.%N.%j.err     # STDERR
#SBATCH --mail-user=ga824@nyu.edu
##SBATCH --dependency=afterok:4708109

##set variables
##REF="/scratch/ga824/hermes/plasmid/GCF_000146045.2_R64_genomic_GAP1_pSG36_wt_HygMX.fna" ###WORKS WITH THIS REFERENCE
REF="/genomics/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GAP1/GCF_000146045.2_R64_genomic_GAP1.fna"

##capture the output of a command line and store it in a variable
ALL_R1=($(ls /scratch/cgsb/gencore/out/Gresham/2019-09-17_HTFWFAFXY/merged/*n01*))

##this is going to assign the variables to file names
FORWARD=${ALL_R1[$SLURM_ARRAY_TASK_ID]}
NAME=${FORWARD:76:-22}

echo $NAME

## change into directory for this specific sample
cd d_$NAME

##INSTEAD of mapping to plasmid, removing those reads, and re-mapping to genome, make the plasmid a contig in the genome

module purge

## map to reference genome with plasmid
## use bwa backtrack, which is optimized for short reads
module purge
module load bwa/intel/0.7.17
## align each read separately
bwa aln -t 3 ${REF} ${NAME}_trimmed_r1.fastq.gz > ${NAME}_aln_r1.sai
bwa aln -t 3 ${REF} ${NAME}_trimmed_r2.fastq.gz > ${NAME}_aln_r2.sai

## create paired end alignment
bwa sampe $REF ${NAME}_aln_r1.sai ${NAME}_aln_r2.sai ${NAME}_trimmed_r1.fastq.gz ${NAME}_trimmed_r2.fastq.gz > ${NAME}_aln_pe.sam

module purge
module load samtools/intel/1.9
##convert to bam
samtools view -bS ${NAME}_aln_pe.sam > ${NAME}_aln_pe.bam

## get stats from alignment
samtools flagstat ${NAME}_aln_pe.bam > ${NAME}_flagstat_bktrk.txt
samtools sort -@ 3 -o ${NAME}_bktk_sort.bam ${NAME}_aln_pe.bam
samtools index ${NAME}_bktk_sort.bam

