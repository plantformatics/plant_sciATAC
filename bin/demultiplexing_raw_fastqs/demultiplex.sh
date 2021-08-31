#!/bin/bash
#source activate cutadaptenv
#bash demultiplex.sh list_demultiplex.config

cat $1 | while read i; do

  # process config file
  set IFS=$"\t"; array=($i); unset IFS
  id=${array[0]}
  read1=${array[1]}${array[2]}
  read2=${array[1]}${array[3]}
  echo "$id $read1 $read2"

  #exit 0

  # find reads with adapter, remove untrimmed
  cutadapt -e 0.1 --no-indels --cores=30 --action=none --discard-untrimmed \
  -g NNNNNNGGCGGTAGGCGTGNNNNNAGATGTGTATAAGAGACAG \
  -G NNNNNNCCAACACCCGTGCGNNNNNAGATGTGTATAAGAGACAG \
  -o ${id}_R1.fastq.gz -p ${id}_R2.fastq.gz  \
  $read1 $read2 &>> ${id}_log.txt

  # select reads matching well barcodes
  cutadapt --cores=30 -e 0.2 --no-indels --discard-untrimmed --action=none --suffix '+{name}' \
  -g well_bc_8row.fa \
  -G well_bc_12column.fa \
  -o ${id}_well_R1.fastq.gz -p ${id}_well_R2.fastq.gz \
  ${id}_R1.fastq.gz ${id}_R2.fastq.gz &>> ${id}_wellbc_log.txt
  rm ${id}_R*.fastq.gz

  # trim well barcode
  cutadapt -u 19 -U 20 \
  -o ./${id}_tmp_R1.fastq.gz -p ./${id}_tmp_R2.fastq.gz \
  ${id}_well_R1.fastq.gz ${id}_well_R2.fastq.gz &>> ${id}_log.txt
  rm ${id}_well_R*.fastq.gz

  # select reads matching cell barcodes
  cutadapt --cores=30 -e 0 --no-indels --discard-untrimmed --action=none --suffix '+{name}' \
  -g cell_bc_24column.fa \
  -G cell_bc_16row.fa \
  -o ${id}_cell_R1.fastq.gz -p ${id}_cell_R2.fastq.gz \
  ${id}_tmp_R1.fastq.gz ${id}_tmp_R2.fastq.gz &>> ${id}_cellbc_log.txt
  rm ${id}_tmp_R*.fastq.gz

  grep "^Sequence:" ${id}_wellbc_log.txt |sed 's/\|;\|://g'| cut -d' ' -f2,9 | awk '{if (NR<=8) {print "scTS"NR" "$0} else {id=NR-8; print "scPCR2."id" "$0} }' > ${id}_report.txt
  grep "^Sequence:" ${id}_cellbc_log.txt |sed 's/\|;\|://g'| cut -d' ' -f2,9 | awk '{if (NR<=24){print "Tn5A"NR" "$0} else {id=NR-24; print "Tn5B"id" "$0} }' >>${id}_report.txt

done
