#!/bin/bash

#SBATCH --partition=highmem_p
#SBATCH --job-name=process_dups_sciATAC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=7-00:00:00
#SBATCH --mem=490g
#SBATCH --output=logs_dedup_sciATAC.%j.log
#SBATCH --error=logs_dedup_sciATAC.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# threads
threads=30

# load modules
ml picard/2.21.6-Java-11

# input files
sci1=Arabidopsis_root_sciATAC_rep1.mq10.bam
sci2=Arabidopsis_root_sciATAC_rep2.mq10.bam

# functions
doCall(){

	# input
	base=$1
	rep=$2
	threads=$3

	# remove duplicates
	echo "removing dups - $base ..."
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		REMOVE_DUPLICATES=true \
		METRICS_FILE=$base.metrics \
		I=$base.mq10.bam \
		O=$base.mq10.BC.rmdup.bam \
		BARCODE_TAG=BC \
		ASSUME_SORT_ORDER=coordinate

	# fix barcodes
	echo "fixing barcodes and removing multi-mapped reads ..."
	perl ./util_scripts/fixBC.pl $base.mq10.BC.rmdup.bam $rep | samtools view -bhS - > $base.mq10.BC.rmdup.mm.bam

	# make Tn5 bed files
	echo "making Tn5 bed files ..."
	perl ./util_scripts/makeTn5bed.pl $base.mq10.BC.rmdup.mm.bam | sort -k1,1 -k2,2n - > $base.tn5.bed

	# remove non-unique Tn5 sites
	uniq $base.tn5.bed > $base.unique.tn5.bed

}
export -f doCall

# run pipeline
doCall Arabidopsis_root_sciATAC_rep1 Arabidopsis_root_sciATAC_rep1 $threads
doCall Arabidopsis_root_sciATAC_rep2 Arabidopsis_root_sciATAC_rep2 $threads
