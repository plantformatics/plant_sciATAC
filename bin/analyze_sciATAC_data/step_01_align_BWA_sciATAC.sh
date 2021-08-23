#!/bin/bash

##----------------------------------------------##
##          SLURM submission properties         ##
##----------------------------------------------##

#SBATCH --partition=highmem_p
#SBATCH --job-name=map_sciATAC_At_root
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=120g
#SBATCH --output=logs_map_sciATAC_At_root.%j.log
#SBATCH --error=logs_map_sciATAC_At_root.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# modules
ml BWA/0.7.17-GCC-8.3.0
ml UMI-tools/1.0.1-foss-2019b-Python-3.7.4

# common variables
fastq=FASTQ
bamout=BAM
ref=TAIR10_reference.fa
threads=24

# rep1
A1=AtR_G019_10x_R1.fastq.gz
A2=AtR_G019_10x_R2.fastq.gz
Ai=AtR_G019_10x_index.fastq.gz
Ai1=AtR_G019_10x_R1.BC.fastq.gz
Ai2=AtR_G019_10x_R2.BC.fastq.gz

# attach barcodes
attachBC(){
	
	# read file
	read1=$1
	
	# index file
	readi=$2

	# output file
	out=$3

	# run
	umi_tools extract --bc-pattern=NNNNNNNNNNNNNNNNNNNNNNNNNNN --stdin=$readi --read2-in=$read1 --stdout=$out --read2-stdout

}
export -f

# attach barcodes
attachBC $fastq/$A1 $fastq/$Ai $fastq/$Ai1
attachBC $fastq/$A2 $fastq/$Ai $fastq/$Ai2

# run BWA
bwa mem -M -t $threads $ref $fastq/$Ai1 $fastq/$Ai2 | samtools view -hbSq 10 | samtools sort - | perl modify_BC_flag.pl - | samtools view -bSh - > Arabidopsis_root_sciATAC_rep1.mq10.bam
