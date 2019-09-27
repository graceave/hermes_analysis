#!/bin/bash                     
#SBATCH --job-name=v4clean_reads
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --array=1-8
#SBATCH --mem=8GB
#SBATCH --mail-type=END
#SBATCH -o reads_clean_relax.%N.%j.out        # STDOUT
#SBATCH -e reads_clean_relax.%N.%j.err     # STDERR
#SBATCH --mail-user=ga824@nyu.edu

##set variables
PREINSERT="^tcataagtagcaagtggcgcataagtatcaaaataagccacttgttgttgttctctg"
BAD="^cccgggggatcc"

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

echo $NAME > cutadapt_new.txt
#report original number of reads
echo "Reads in library, before cleaning:" >> cutadapt_new.txt
echo $(zcat ${FORWARD}|wc -l)/4|bc >> cutadapt_new.txt


 ## 2 successive cleaning steps
    ## - take out UMIs, associate with read name
    ## - trim the plasmidic sequence before the insertion point, remove the untrimmed reads 

 module purge

## extract UMIs
module load umi_tools/0.5.5
umi_tools extract --stdin=${FORWARD} --read2-in=${REVERSE} --bc-pattern=NNNNNNN --log=umi_extract.log --stdout ${NAME}_umi_extracted_r1.fastq.gz --read2-out=${NAME}_umi_extracted_r2.fastq.gz

module purge

## trim expected plasmid off
## only keep reads where there was the preinsert sequence
## increase maximum error rate to 30%
## require Minimum overlap 40 bp
module load cutadapt/1.16
cutadapt -g ${PREINSERT} --discard-untrimmed -e 0.2 -O 40 -o ${NAME}_trimmed_r1.fastq.gz -p ${NAME}_trimmed_r2.fastq.gz ${NAME}_umi_extracted_r1.fastq.gz ${NAME}_umi_extracted_r2.fastq.gz >> cutadapt_new.txt

## remove any reads that have the plasmid sequence continuing after the preinsert
cutadapt -g ${BAD} --discard-trimmed -o ${NAME}_clean_r1.fastq.gz -p ${NAME}_clean_r2.fastq.gz ${NAME}_trimmed_r1.fastq.gz ${NAME}_trimmed_r2.fastq.gz >> cutadapt_new.txt

##HERE report number of reads after searching for the pre-insert sequence
echo "Reads in library (forward), after trimming:" >> cutadapt_new.txt
echo $(zcat ${NAME}_clean_r1.fastq.gz|wc -l)/4|bc >> cutadapt_new.txt
echo "Reads in library (reverse), after trimming:" >> cutadapt_new.txt
echo $(zcat ${NAME}_clean_r2.fastq.gz|wc -l)/4|bc >> cutadapt_new.txt


