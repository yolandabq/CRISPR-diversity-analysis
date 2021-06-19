#!/bin/bash

# Autora: Yolanda Benítez 
# Template script to build the fastq from the lab data and obtain the Bam files (to be analyzed with CrispRVariants).

fastq_list=$(ls *[0-9].fastq | cut -f 1 -d ".")
fastq_list_uq=$(echo $fastq_list | tr ' ' '\n' | sort -nu)

for f in $fastq_list_uq 
do

  #echo $f
  cat $f*.fastq > ./fastq_merged/$f.fastq
  	
done

# Look at the qualities of the fastq

for f in $fastq_list_uq 
do
  	
  bwa mem /mnt/c/Users/yolib/Documents/mm9_index_bwa/mm9 ./fastq_merged/$f.fastq > ./fastq_merged/$f.sam
  samtools view -S -b ./fastq_merged/$f.sam > ./fastq_merged/$f.bam 


done

