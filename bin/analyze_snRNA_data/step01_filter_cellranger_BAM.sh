#!/bin/bash

# submission properties --------------------------------------
#SBATCH --partition=highmem_p
#SBATCH --job-name=clean_At_root_10xsnRNA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6-23:00:00
#SBATCH --mem=490g
#SBATCH --output=logs_clean_At_root_10xsnRNA.%j.log
#SBATCH --error=logs_clean_At_root_10xsnRNA.%j.err

# pre-launch settings ----------------------------------------
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# parameters
ncpu=1
mem=490g

# modules
ml picard/2.21.6-Java-11

# directories
indir=$PWD

# input
rep1=Atroot_snRNAseq_rep1
rep2=Atroot_snRNAseq_rep2
rep3=Atroot_snRNAseq_rep3
rep4=Atroot_snRNAseq_rep4
rep5=Atroot_snRNAseq_rep5

# functions -------------------------------------------------
dedup(){

        # variables
        i=$1
        o=$2
        m=$3
        mem=$4

	# state parameters
	echo "input = $i | output = $o | metrics = $m | memory = $mem"

        # run
        java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        REMOVE_DUPLICATES=true \
        METRICS_FILE=$m \
        I=$i \
        O=$o \
        BARCODE_TAG=UM \
        ASSUME_SORT_ORDER=coordinate

}
export -f dedup

# run -------------------------------------------------------
perl ./utils/combine_CB_UB.pl $rep1.raw.bam 1 | samtools view -bS - > $rep1.v1.bam
perl ./utils/combine_CB_UB.pl $rep2.raw.bam 1 | samtools view -bS - > $rep2.v1.bam
perl ./utils/combine_CB_UB.pl $rep3.raw.bam 1 | samtools view -bS - > $rep3.v1.bam
perl ./utils/combine_CB_UB.pl $rep4.raw.bam 1 | samtools view -bS - > $rep4.v1.bam
perl ./utils/combine_CB_UB.pl $rep5.raw.bam 1 | samtools view -bS - > $rep5.v1.bam
dedup $rep1.v1.bam $rep1.v2.bam $rep1.metrics $mem 
dedup $rep2.v1.bam $rep2.v2.bam $rep2.metrics $mem 
dedup $rep3.v1.bam $rep3.v2.bam $rep3.metrics $mem 
dedup $rep4.v1.bam $rep4.v2.bam $rep4.metrics $mem 
dedup $rep5.v1.bam $rep5.v2.bam $rep5.metrics $mem
perl ./utils/convertBAM2sparse.pl $rep1.v2.bam Atroot_sn_rep1 > $rep1.sparse
perl ./utils/convertBAM2sparse.pl $rep2.v2.bam Atroot_sn_rep2 > $rep2.sparse
perl ./utils/convertBAM2sparse.pl $rep3.v2.bam Atroot_sn_rep3 > $rep3.sparse
perl ./utils/convertBAM2sparse.pl $rep4.v2.bam Atroot_sn_rep4 > $rep4.sparse
perl ./utils/convertBAM2sparse.pl $rep5.v2.bam Atroot_sn_rep5 > $rep5.sparse
perl ./utils/countData.pl $rep1.v2.bam Atroot_sn_rep1 > $rep1.metadata.txt
perl ./utils/countData.pl $rep2.v2.bam Atroot_sn_rep2 > $rep2.metadata.txt
perl ./utils/countData.pl $rep3.v2.bam Atroot_sn_rep3 > $rep3.metadata.txt
perl ./utils/countData.pl $rep4.v2.bam Atroot_sn_rep4 > $rep4.metadata.txt
perl ./utils/countData.pl $rep5.v2.bam Atroot_sn_rep5 > $rep5.metadata.txt

# merge
cat *.sparse > Atroot_snRNAseq_merged.sparse
cat *.metadata.txt > Atroot_snRNAseq_merged.metadata.raw.txt
head -n 1 Atroot_snRNAseq_merged.metadata.raw.txt > header
grep -v 'cellID' Atroot_snRNAseq_merged.metadata.raw.txt > temp.meta
cat header temp.meta > Atroot_snRNAseq_merged.metadata.raw.txt
rm header temp.meta
