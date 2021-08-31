#!/bin/bash
#merge demultiplex read R1 R2 from each plates (still have the cell bc at the 5' end) to fastq of 10X format: R1, index read and R2
#a 5nt bc is added to distinguish reads from each plate
#input = name + 5nt barcode + file1 + file2
#bash combine_10x.sh list_combine.config

cat $1 | while read i; do
set IFS=$"\t"; array=($i); unset IFS
id=${array[0]}
bc=${array[1]}
read1=${array[2]}
read2=${array[3]}

paste <(zcat ${read1}) <(zcat ${read2}) | awk -v awk_id=$id -v awk_bc=$bc '{
q="JJJJJFFFFFFFFFFFFFFFFFFFFFF"
f1=awk_id"_10x_R1.fastq.gz"
f2=awk_id"_10x_R2.fastq.gz"
f3=awk_id"_10x_index.fastq.gz"
if (NR%4==1) {
	split($2,index1,/+/)
	split($4,index2,/+/)
	barcode=index1[3] index1[4] index2[3] index2[4]
	print $1" "$2"+"index2[3]"+"index2[4] | "gzip >>" f1
	sub(/2:N:0/,"3:N:0",$4)	#change the read 2 header label to read3, while index read become read2
	print $3" "$4 | "gzip >>" f2
	sub(/1:N:0/,"2:N:0",$2)
	print $1" "$2 | "gzip >>" f3;}
else if (NR%4==2) {
        print substr($1,24) | "gzip >>" f1
        print substr($2,24) | "gzip >>" f2
        print awk_bc barcode | "gzip >>" f3;}
else if (NR%4==3) {
        print $1 | "gzip >>" f1
        print $2 | "gzip >>" f2
        print $1 | "gzip >>" f3;}
else  {
        print q | "gzip >>" f3
        print substr($1,24) | "gzip >>" f1
        print substr($2,24) | "gzip >>" f2;}

}'

done
