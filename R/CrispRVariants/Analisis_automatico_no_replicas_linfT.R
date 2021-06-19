# Script to analyze T lymphocyte data with CrispRVariants.
# Analysis of the guides with less than 3 replicates These go to the "result_no_repl" and "result_no_replic_ins" folders.


library(CrispRVariants)
library("rtracklayer")
library("GenomicFeatures")
library("gdata")
library(tidyverse)
library("RColorBrewer")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library("GenomicAlignments")
library("Rsamtools")
library(dplyr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


setwd("C:/Users/yolib/Documents/TFM/Linf_T/Datos")

####################################

#files_info <- read.csv("replicas.txt", sep = "\t", header = TRUE)

all_files_info <- read.csv("all_target_replicas.txt", sep = "\t", header = TRUE)
guide_data <- read.csv("guideData.csv", sep = "\t")

unique_guides <- data.frame(guide_name = character(),
                            bam_files = character())


# I organize a data frame with the name of the guides and their corresponding bams.
# The guides separated by ",", the bams of each guide, separated by ";".

for (gene in unique(guide_data$genename)){
  
  guides <- guide_data$guide[guide_data$genename == gene]
  
  bam_str = ''
  
  for (guide in guides){
    #bam_str <- str_c(bam_str, ',', all_files_info$Bam_names[all_files_info$Guia_ID == guide])
    if (!identical(all_files_info$Bam_names[all_files_info$Guia_ID ==  guide], character(0))){
      bam_str <- str_c(bam_str, ',', substr(all_files_info$Bam_names[all_files_info$Guia_ID == guide], 2, nchar(all_files_info$Bam_names[all_files_info$Guia_ID == guide])))
    }
    else {
      guides <- guides[-c(which(guides == guide))]
      
    }
  }
  
  # I will save the ones that have more than 3 bams
  
  bam_str <- substr(bam_str, 2, nchar(bam_str))
  guides_name <-paste(guides, collapse=',' ) 
  # the "," separates the bam from different guides, and the ";" separate the bams within the same guide
  
  
  #if (bam_str != ''){
  if (length(strsplit(str_replace_all(bam_str, ";", ","), ",")[[1]]) < 3){
    unique_guides <- rbind(unique_guides, c(guides_name, bam_str))
  }
  
  
}

colnames(unique_guides) <- c("Guides_names", "Bam_files")

#View(unique_guides[0:15,])

####################################


# This part is importante to obtain the ID and the list of bam.
# The bams of the same guide are separated by ";", and those belonging to different guides, by ",".

###

n_replic = nrow(unique_guides) # number of replicas (of different guides)

ref_data = data.frame(guide_id = character(),
                      sequence = character())
#guide_index = 0

for (replic in c(1:n_replic)){
  #for (replic in c(1:5)){
  
  #guide_id_file <- files_info[replic,]$Guia_ID
  
  list_guide_id_file <- strsplit(str_replace_all(unique_guides[replic,]$Guides_names, ",", ","), ",")[[1]]
  #guide_id_name <- strsplit(guide_id_file, '-')[[1]][2]
  #bams_id <- unique_guides[replic,]$Bam_files
  bams_id_pergRNA <- strsplit(str_replace_all(unique_guides[replic,]$Bam_files, ",", ","), ",")[[1]]
  #bams_id <- substr(bams_id, 2, nchar(bams_id))
  
  g = 1
  #guide_index = guide_index + 1
  
  for (guide_id_name in list_guide_id_file){
    
    setwd("C:/Users/yolib/Documents/TFM/Linf_T/Datos/Datos")
    
    #print(guide_id_name)
    
    gRNA_files_name = list.files(pattern = str_c("-", guide_id_name))
    
    for (guide_file_name in gRNA_files_name){
      
      #print(guide_file_name)
      
      guide_path <- str_c("C:/Users/yolib/Documents/TFM/Linf_T/Datos/Datos/", guide_file_name)
      setwd(guide_path)
      
      
      
      list_bam_full <- strsplit(str_replace_all(bams_id_pergRNA[g], ";", ","), ",")[[1]]
      #print(list_bam_full)
      
      g = g + 1
      names <- str_replace(list_bam_full, ".bam", "")
      
      reference <-  guide_data$reference[guide_data$names == guide_id_name]
      start <- as.numeric(guide_data$starts[guide_data$names == guide_id_name])
      end <- as.numeric(guide_data$ends[guide_data$names == guide_id_name])
      chr <- guide_data$seqnames[guide_data$names == guide_id_name]
      strand <- guide_data$strands[guide_data$names == guide_id_name]
      gene_name <- guide_data$genename[guide_data$names == guide_id_name]
      
      ref <- Biostrings::DNAString(reference) # voy a utilizar la mutada
      gd <- GenomicRanges::GRanges(chr, IRanges::IRanges(start + 1, end), strand = strand, genome = "hg38")
      
      
      ## -------------------------- DECOMMENT ALL OF THIS TO EXECUTE IT ----------------------------
      crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref, 
                                  target.loc = 30,  names = names) #collapse.pairs = TRUE, no lo pongo porque ellos fusionaron las reads con FLASH
      ## ---------------------------------------------------------------------------------------------
      ref_seq <- crispr_set$ref
      
      
      ref_data <- rbind(ref_data, c(as.character(replic), as.character(ref_seq)))
      
      
      # If we only want the alleles with a frequency > X%: 
      
      #common_alleles <- as.data.frame(variantCounts(crispr_set, result = "proportions", min.freq = 1),
      #row.names = rownames(variantCounts(crispr_set, result = "proportions", min.freq = 1)))
      
      ## -------------------------- DECOMMENT ALL OF THIS TO EXECUTE IT----------------------------
      # If we want all the alleles: 
      common_alleles <- as.data.frame(variantCounts(crispr_set, result = "proportions"),
                                      row.names = rownames(variantCounts(crispr_set, result = "proportions")))
      
      common_alleles['Alleles'] = rownames(common_alleles) # para guardar el nombre de la columna 
      
      
      file_path <- str_c("C:/Users/yolib/Documents/TFM/Linf_T/Datos/results_no_replic/", replic + 526, "_", guide_id_name, ".csv") # para guardar el csv
      
      write.csv(common_alleles, file_path, row.names = FALSE)
      
      insertions <- as.data.frame(crispr_set$insertion_sites[, c(1,2,3,4,6,7)])
      
      insertions["sample_name"] =  names[insertions$sample]
      
      insertions["proportion"] = unlist(lapply(insertions$count, function(x) {100*x/sum(insertions$count)}))
      
      file_path <- str_c("C:/Users/yolib/Documents/TFM/Linf_T/Datos/results_no_replic_ins/", replic + 526, "_", guide_id_name, "_ins.csv") # para guardar el csv
      
      write.csv(insertions, file_path, row.names = FALSE)
      # ----------------------------------------------------------------------------------------------
    }
    
    ref_data <- rbind(ref_data, c(as.character(replic + 526), as.character(reference), as.character(strand)))
  }
  
  colnames(ref_data) = c("Guide_ID", "Sequence", "Strand")
  file_path <- str_c("C:/Users/yolib/Documents/TFM/Linf_T/Datos/guides_sequence_2_no_replic.csv") # para guardar el csv
  write.csv(ref_data, file_path, row.names = FALSE)
}

'#
  list_bam_full <- strsplit(str_replace_all(bams_id, ";", ","), ",")[[1]]
  names <- str_replace(list_bam_full, ".bam", "")
  
  
  reference <-  guide_data$reference[guide_data$names == guide_id_name]
  start <- as.numeric(guide_data$starts[guide_data$names == guide_id_name])
  end <- as.numeric(guide_data$ends[guide_data$names == guide_id_name])
  chr <- guide_data$seqnames[guide_data$names == guide_id_name]
  strand <- guide_data$strands[guide_data$names == guide_id_name]
  gene_name <- guide_data$genename[guide_data$names == guide_id_name]
  
  guide_path <- str_c("C:/Users/yolib/Documents/TFM/Linf_T/Datos/Datos/", guide_id_file)
  setwd(guide_path)
  #list_bam_full <- list.files(pattern = "*[0-9].bam$")
  
  ref <- Biostrings::DNAString(reference) # voy a utilizar la mutada
  gd <- GenomicRanges::GRanges(chr, IRanges::IRanges(start + 1, end), strand = strand, genome = "hg38")
  
  
  
  crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref, 
                              target.loc = 30,  names = names) #collapse.pairs = TRUE, no lo pongo porque ellos fusionaron las reads con FLASH
  
  #p <- plotVariants(crispr_set, top.plot = 1,  
  #                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency")) # ahora el header es la eficiencia 
  
  #### AHORA VOY A GUARDAR LOS DATOS EN ARCHIVOS PARA PODER ANALIZARLOS CON PYTHON
  
  ref_seq <- crispr_set$ref
  
  #common_alleles <- as.data.frame(variantCounts(crispr_set, result = "proportions"),
                                  #row.names = rownames(variantCounts(crispr_set, result = "proportions")))
  
  #common_alleles <- as.data.frame(variantCounts(crispr_set, result = "proportions", min.freq = 1),
                                 #row.names = rownames(variantCounts(crispr_set, result = "proportions", min.freq = 1)))
  
  common_alleles <- as.data.frame(variantCounts(crispr_set,  min.freq = 1),
                                  row.names = rownames(variantCounts(crispr_set, min.freq = 1)))
  
  
  # ordeno los alelos de más a menos frecuente (ya que crisprvariants pone primero los SNVs, aunque tengan menos porcentaje)
  
  common_alleles["Alleles"] = rownames(common_alleles)
  file_path <- str_c("C:/Users/yolib/Documents/TFM/Linf_T/Datos/results_min_freq1/", guide_id_name, "_counts.csv")
  
  write.csv(common_alleles, file_path, row.names = FALSE)
  
  
}
 
'
 

