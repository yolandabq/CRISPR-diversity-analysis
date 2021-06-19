#!/bin/bash

# Script to preprocess SRA data from T lymphocytes (Ryan T. Leenay et al.)

sra_list=$(cut -f 2 -d "," /mnt/c/Users/yolib/Documents/sra_1000_paired.csv | tail -n +2 | sed -e 's/^.//' -e 's/.$//' | head -n 2)

echo $sra_list


for sra_i in $sra_list
do
 echo $sra_i # we create the folders
 if [ ! -d "/mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i" ]; then 
	mkdir /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i
 fi 
 
 prefetch $sra_i -O /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i  # we download the sra file
 fastq-dump -I --split-3 -O /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/$sra_i.sra # we convert it in fastq
 
 # we remove the adapters and check the quality. /mnt/c/Users/yolib/OneDrive/Documentos/Software/Trimmomatic-0.39
 
 
  if [ -e /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_2.fastq ]; then # if it is paired end
   
   /mnt/c/Users/yolib/OneDrive/Documentos/Software/FLASH-1.2.11-Linux-x86_64/flash /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_1.fastq /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_2.fastq -d /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/ -o $sra_i -q
   
   ### -----------------------------------------WITH FLASH WE JOIN THE TWO PAIRED END ----------------------------------------------------------
   #java -jar /mnt/c/Users/yolib/OneDrive/Documentos/Software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_1.fastq /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_2.fastq /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_1_paired.fastq /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_1_unpaired.fastq /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_2_paired.fastq /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_2_unpaired.fastq ILLUMINACLIP:/mnt/c/Users/yolib/OneDrive/Documentos/Software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20  ## If we do trimmomatic before flash 
   
   java -jar /mnt/c/Users/yolib/OneDrive/Documentos/Software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i".extendedFrags.fastq /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_trimmed.fastq SLIDINGWINDOW:4:20 ILLUMINACLIP:/mnt/c/Users/yolib/OneDrive/Documentos/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 ## If we do flash before trimmomatic 
   ## java -jar /mnt/c/Users/yolib/OneDrive/Documentos/Software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 /mnt/c/Users/yolib/Documents/TFM/Linf_T/SRR7752829/SRR7752829.extendedFrags.fastq /mnt/c/Users/yolib/Documents/TFM/Linf_T/SRR7752829/SRR7752829_trimmed.fastq SLIDINGWINDOW:4:20 ILLUMINACLIP:/mnt/c/Users/yolib/OneDrive/Documentos/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 
   
  else # if it is single end
   #echo "single"
   java -jar /mnt/c/Users/yolib/OneDrive/Documentos/Software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i".fastq /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_trimmed.fastq SLIDINGWINDOW:4:20 ILLUMINACLIP:/mnt/c/Users/yolib/OneDrive/Documentos/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10
  
  fi
 
 # alignment with BWA 
 
 bwa mem /mnt/c/Users/yolib/Documents/hg38_index_bwa/hg38.fa.gz /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i"_trimmed.fastq > /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i".sam
 ## bwa mem /mnt/c/Users/yolib/Documents/hg38_index_bwa/hg38.fa.gz /mnt/c/Users/yolib/Documents/TFM/Linf_T/SRR7752829/SRR7752829_trimmed.fastq > /mnt/c/Users/yolib/Documents/TFM/Linf_T/SRR7752829/SRR7752829.sam

 
 samtools view -S -b /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i".sam > /mnt/c/Users/yolib/Documents/TFM/Linf_T/$sra_i/"$sra_i".bam 
##samtools view -S -b /mnt/c/Users/yolib/Documents/TFM/Linf_T/SRR7752829/SRR7752829.sam > /mnt/c/Users/yolib/Documents/TFM/Linf_T/SRR7752829/SRR7752829.bam 

