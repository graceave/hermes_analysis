#!/bin/bash                     
#SBATCH --job-name=v5clean_reads
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --array=0
#SBATCH --mem=8GB
#SBATCH --mail-type=END
#SBATCH -o reads_clean_relax.%N.%j.out        # STDOUT
#SBATCH -e reads_clean_relax.%N.%j.err     # STDERR
#SBATCH --mail-user=ga824@nyu.edu

##set variables
PREINSERT="^tcataagtagcaagtggcgcataagtatcaaaataagccacttgttgttgttctctg"
BAD="^ggatcccccggg"

##capture the output of a command line and store it in a variable
ALL_R1=($(ls /scratch/cgsb/gencore/out/Gresham/2019-09-17_HTFWFAFXY/merged/*n01*))
ALL_R2=($(ls /scratch/cgsb/gencore/out/Gresham/2019-09-17_HTFWFAFXY/merged/*n02*))

##this is going to assign the variables to file names
FORWARD=${ALL_R1[$SLURM_ARRAY_TASK_ID]}
REVERSE=${ALL_R2[$SLURM_ARRAY_TASK_ID]}
NAME=${FORWARD:76:-22}

##change to directory for this specific sample
mkdir d_$NAME
cd d_$NAME



 ## 3 successive cleaning steps
    ## - take out UMIs, associate with read name
    ## - trim the plasmidic sequence before the insertion point, remove the untrimmed reads
    ## - trim the plasmidic sequence after the insertion point, remove the trimmed reads 

 module purge

## extract UMIs
module load umi_tools/0.5.5
umi_tools extract --stdin=${FORWARD} --read2-in=${REVERSE} --bc-pattern=NNNNNNN --log=umi_extract.log --stdout ${NAME}_umi_extracted_r1.fastq.gz --read2-out=${NAME}_umi_extracted_r2.fastq.gz

module purge

## trim expected plasmid off
## only keep reads where there was the preinsert sequence
## increase maximum error rate to 30%
## require Minimum overlap 40 bp
## require minimum read length of 10 bases
module load cutadapt/1.16
cutadapt -g ${PREINSERT} --discard-untrimmed -e 0.3 -O 40 -m 10 -o ${NAME}_trimmed_r1.fastq.gz -p ${NAME}_trimmed_r2.fastq.gz ${NAME}_umi_extracted_r1.fastq.gz ${NAME}_umi_extracted_r2.fastq.gz > cutadapt_new.txt


## remove any reads that have the plasmid sequence continuing after the preinsert
cutadapt -g ${BAD} --discard-trimmed -e 0 -o ${NAME}_clean_r1.fastq.gz -p ${NAME}_clean_r2.fastq.gz ${NAME}_trimmed_r1.fastq.gz ${NAME}_trimmed_r2.fastq.gz >> cutadapt_new.txt

##HERE report number of reads at each stage
echo $NAME > ${NAME}_cleaning_counts.txt
echo "Reads in library (f), before cleaning:" >> ${NAME}_cleaning_counts.txt
echo $(zcat ${FORWARD}|wc -l)/4|bc >> ${NAME}_cleaning_counts.txt
echo "Reads in library, with hermes "preinsert":" >> ${NAME}_cleaning_counts.txt
echo $(zcat ${NAME}_trimmed_r1.fastq.gz|wc -l)/4|bc >> ${NAME}_cleaning_counts.txt
echo "Reads in library (forward), after trimming:" >> ${NAME}_cleaning_counts.txt
echo $(zcat ${NAME}_clean_r1.fastq.gz|wc -l)/4|bc >> ${NAME}_cleaning_counts.txt
echo "Reads in library (reverse), after trimming:" >> ${NAME}_cleaning_counts.txt
echo $(zcat ${NAME}_clean_r2.fastq.gz|wc -l)/4|bc >> ${NAME}_cleaning_counts.txt


