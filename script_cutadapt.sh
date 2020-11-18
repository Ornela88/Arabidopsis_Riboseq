#!/bin/bash
#add adapter sequence you want to trim
adapter_sequence=TGGAATTCTCGGGTGCCAAGG
home=/usr/data/bgfs1/maloku/

#Add the directory containing the raw fastq reads (raw reads should be gzipped)
raw_data=/usr/data/bgfs1/maloku/fastq-raw_data/

####################remove extension tag of files and store sample names in a file for easy access

for fq in ${raw_data}/*fastq.gz; do
	    sample=$(echo $fq | sed 's/.*\///' | sed 's/\.fastq\.gz//')
	            echo $sample >> ${home}/sample_name
	    done


################################### #Trimming using cutadapt

for file in $(cat ${home}/sample_name)
do
	        echo $file "trimming"
		        (cutadapt --cores=8 --minimum-length 25 -a $adapter_sequence -o ${file}_trimmed.fastq.gz ${raw_data}/${file}.fastq.gz) 1> ${file}_trimmed_log
				       
done


