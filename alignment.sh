#!/bin/bash

home=/usr/data/bgfs1/maloku/
#Add the directory containing the raw fastq reads (raw reads should be gzipped)
trimm_data=/usr/data/bgfs1/maloku/trimm/
ribo_bam=/usr/data/bgfs1/maloku/ribo_bam/


######################## Unzip *trimmed.fastq.gz ###########

cd trimm_data/
gzip -d *.fastq.gz


####################remove extension tag of files and store sample names in a file for easy access

for fq in ${trimm_data}/*trimmed.fastq; 
do
		   
 sample=$(echo $fq | sed 's/.*\///' | sed 's/\_trimmed\.fastq//')
	echo $sample >> ${home}/sample_name_trimmed
done
															       
###################### HISAT2 for  alignment of reads on transcriptome

## Create an index file

hisat2-build -f /tair10/Arabidopsis_thaliana.TAIR10.cdna.all.fa

mkdir index
mv Athaliana* index/

##### Align reads and create SAM files

for file in $(cat ${home}/sample_name_trimmed);
do	                              
  echo $file "Align reads"

  (hisat2 -p 8 --dta -x index/Athaliana -q ${trimm_data}/${file}_trimmed.fastq -S ${ribo_bam}/${file}.sam)
done


##### convert output SAM files to sorted BAM files

for file in $(cat ${home}/sample_name_trimmed);
do
    echo $file "Create BAM files"

samtools sort -@ 4 -o  ${ribo_bam}/${file}.bam  ${ribo_bam}/${file}.sam
done

#### Create the .BAI files from BAM files in the same directory

cd ribo_bam/

for i in *.bam

do

	echo "Indexing: "$i        

	samtools index $i $i".bai"

done

