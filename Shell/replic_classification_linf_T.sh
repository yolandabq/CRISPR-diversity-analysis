#!/bin/bash
# This script is to obtain the name of the bam that belong to each file (with more than 1000 reads). It create a txt called replicas with the names of the corresponding guides.

echo "Guia_ID	Bam_names"  > /mnt/c/Users/yolib/Documents/TFM/Linf_T/Datos/all_target_replicas.txt 

cd /mnt/c/Users/yolib/Documents/TFM/Linf_T/Datos/Datos

carpetas=$(ls)

#echo $carpetas
for file in $carpetas # We enter into the folders with the data one by one 
do 

  guide_name=$(echo $file | cut -d '-' -f2)
  
  cd $file # we enter into the folder
  
  bam_files=$(ls | grep '.bam$') # we obtain the name of the .bam
  cadena_bam=$(echo '') # we initiate the variable "cadena_bam"
  for bam_file in $bam_files  # We go through the bam files that are in each folder
	  do
	  	n_alignm=$(samtools view $bam_file | wc -l) # I see how many alignments each bam has. They do not have header
	  	if (("$n_alignm" >= "1000")) # if it has more than 1000 alignments we save it in string_bam. The name of each file is separated by ";"
  		then
  			#echo $bam_file
  			cadena_bam=$(echo $cadena_bam";"$bam_file)
  			
  		fi
  		
	  done
   
   if ! [ -z $cadena_bam ] # If the chain is not empty
   then 
   	echo "$guide_name	$cadena_bam"  >> /mnt/c/Users/yolib/Documents/TFM/Linf_T/Datos/all_target_replicas.txt 
   	#echo $guide_name
   fi
   cd ..
   
done
