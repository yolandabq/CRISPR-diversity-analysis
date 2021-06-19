
# Script to compare the analysis results with CrispRVariants when I generate the files
# .bam with both fastqs, when I do it only with one of the fastqs of the pair and when I use the files
# .bam already processed by Ryan T. Leenay et al.


library(CrispRVariants)
library("rtracklayer")
library("GenomicFeatures")
library("gdata")
library(tidyverse)
library("RColorBrewer")

###############################################################################################

# Bam obtained manually from the raw data 

setwd("C:/Users/yolib/Documents/TFM/Linf_T/SRR7762567")

list_bam_full <- list.files(pattern = '*[0-9].bam$')


ref <- Biostrings::DNAString("GAAGCTTCTGCTCCTATCCCTCACGATGGAAGTAGGTTTATACTTTGAGTTCATTCTCTG") 
gd <- GenomicRanges::GRanges("chr8", IRanges::IRanges(58570935, 58570994), strand = '+', genome = "hg38")

names <- str_replace(list_bam_full, ".bam", "")


crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 30, collapse.pairs = TRUE, names = names)




#crispr_set_opposite <- readsToTarget(list_bam_full, target = gd, reference = ref,
#target.loc = 17, collapse.pairs = TRUE, orientation = "opposite")

p <- plotVariants(crispr_set, top.plot = 1, plot.title =  "F1", 
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency")) 

# Analysis with only one of the fastq of the pair. 

setwd("C:/Users/yolib/Documents/TFM/Linf_T/SRR7762567/single_pairing/")

list_bam_full <- list.files(pattern = '*[0-9].bam$')


ref <- Biostrings::DNAString("GAAGCTTCTGCTCCTATCCCTCACGATGGAAGTAGGTTTATACTTTGAGTTCATTCTCTG") # voy a utilizar la mutada
gd <- GenomicRanges::GRanges("chr8", IRanges::IRanges(58570935, 58570994), strand = '+', genome = "hg38")

names <- str_replace(list_bam_full, ".bam", "")


crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 30, collapse.pairs = TRUE, names = names)

p <- plotVariants(crispr_set, top.plot = 1, plot.title =  "F1", 
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency")) # ahora el header es la eficiencia 




## 640-1_F04. Bam preprocessed by the author. 

setwd("C:/Users/yolib/Documents/TFM/Linf_T/Datos")

files_info <- read.csv("all_target_replicas.txt", sep = "\t", header = TRUE)
guide_data <- read.csv("guideData.csv", sep = "\t")




reference <-  guide_data$reference[guide_data$names == "1_F04"]
start <- guide_data$starts[guide_data$names == "1_F04"]
end <- guide_data$ends[guide_data$names == "1_F04"]
chr <- guide_data$seqnames[guide_data$names == "1_F04"]
strand <- guide_data$strands[guide_data$names == "1_F04"]


setwd("C:/Users/yolib/Documents/TFM/Linf_T/Datos/Datos/640-1_F04")
list_bam_full <- list.files(pattern = '*[0-9].bam$')

bam <- "RL384-00022_K08.bam"
name <- "RL384-00022_K08"

ref <- Biostrings::DNAString(reference)
gd <- GenomicRanges::GRanges(chr, IRanges::IRanges(start + 1, end), strand = strand, genome = "hg38")

names <- str_replace(list_bam_full, ".bam", "")

crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 30, collapse.pairs = TRUE, names = names)

p <- plotVariants(crispr_set, top.plot = 1,  
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency")) 


####################################################################################################################


setwd("C:/Users/yolib/Documents/TFM/Linf_T/SRR7762568")

list_bam_full <- list.files(pattern = '*[0-9].bam$')


ref <- Biostrings::DNAString("AACTGAAGAAAAGATGTCCCTGTACGATGACCTAGGAGTGGAGACCAGTGACTCAAAAAC") 
gd <- GenomicRanges::GRanges("chr10", IRanges::IRanges(6097053, 6097112), strand = '+', genome = "hg38")

names <- str_replace(list_bam_full, ".bam", "")


crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 30, collapse.pairs = TRUE, names = names)

#crispr_set_opposite <- readsToTarget(list_bam_full, target = gd, reference = ref,
#target.loc = 17, collapse.pairs = TRUE, orientation = "opposite")

p <- plotVariants(crispr_set, top.plot = 1, plot.title =  "F1", 
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency"))  



# Este es con solo una de las lecturas del pair end. 

setwd("C:/Users/yolib/Documents/TFM/Linf_T/SRR7762568/single_pairing/")

list_bam_full <- list.files(pattern = '*[0-9].bam$')


ref <- Biostrings::DNAString("AACTGAAGAAAAGATGTCCCTGTACGATGACCTAGGAGTGGAGACCAGTGACTCAAAAAC") 
gd <- GenomicRanges::GRanges("chr10", IRanges::IRanges(6097053, 6097112), strand = '+', genome = "hg38")

names <- str_replace(list_bam_full, ".bam", "")


crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 30, collapse.pairs = TRUE, names = names)

p <- plotVariants(crispr_set, top.plot = 1, plot.title =  "F1", 
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency")) 






### 597-1_B09
setwd("C:/Users/yolib/Documents/TFM/Linf_T/Datos")

files_info <- read.csv("all_target_replicas.txt", sep = "\t", header = TRUE)
guide_data <- read.csv("guideData.csv", sep = "\t")

reference <-  guide_data$reference[guide_data$names == "1_B09"]
start <- guide_data$starts[guide_data$names == "1_B09"]
end <- guide_data$ends[guide_data$names == "1_B09"]
chr <- guide_data$seqnames[guide_data$names == "1_B09"]
strand <- guide_data$strands[guide_data$names == "1_B09"]


setwd("C:/Users/yolib/Documents/TFM/Linf_T/Datos/Datos/597-1_B09")
list_bam_full <- list.files(pattern = '*[0-9].bam$')


bam <- "RL384-00022_C18.bam"
name <- "RL384-00022_C18"

ref <- Biostrings::DNAString(reference) 
gd <- GenomicRanges::GRanges(chr, IRanges::IRanges(start + 1, end), strand = strand, genome = "hg38")

names <- str_replace(list_bam_full, ".bam", "")

crispr_set <- readsToTarget(bam, target = gd, reference = ref,
                            target.loc = 30, collapse.pairs = TRUE, names = name)

p <- plotVariants(crispr_set, top.plot = 1,  
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency"))  



############################################################################33


setwd("C:/Users/yolib/Documents/TFM/Linf_T/SRR7762569")

list_bam_full <- list.files(pattern = '*[0-9].bam$')


ref <- Biostrings::DNAString("CTCCTTTTATCTGTTCTTTGGTACCAGAAGATACGGTTGTTACTGCTGTAGACTGTATCA") 
gd <- GenomicRanges::GRanges("chr4", IRanges::IRanges(106328139, 106328198), strand = '-', genome = "hg38")

names <- str_replace(list_bam_full, ".bam", "")




crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 30, collapse.pairs = TRUE, names = names)

#crispr_set_opposite <- readsToTarget(list_bam_full, target = gd, reference = ref,
#target.loc = 17, collapse.pairs = TRUE, orientation = "opposite")

p <- plotVariants(crispr_set, top.plot = 1, plot.title =  "F1", 
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency")) 




# Este es con solo una de las lecturas del pair end. 

setwd("C:/Users/yolib/Documents/TFM/Linf_T/SRR7762569/single_pairing/")

list_bam_full <- list.files(pattern = '*[0-9].bam$')


ref <- Biostrings::DNAString("CTCCTTTTATCTGTTCTTTGGTACCAGAAGATACGGTTGTTACTGCTGTAGACTGTATCA") 
gd <- GenomicRanges::GRanges("chr4", IRanges::IRanges(106328139, 106328198), strand = '-', genome = "hg38")

names <- str_replace(list_bam_full, ".bam", "")


crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 30,  names = names)

p <- plotVariants(crispr_set, top.plot = 1, plot.title =  "F1", 
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency"))





##612-1_C12


files_info <- read.csv("all_target_replicas.txt", sep = "\t", header = TRUE)
guide_data <- read.csv("guideData.csv", sep = "\t")

reference <-  guide_data$reference[guide_data$names == "1_C12"]
start <- guide_data$starts[guide_data$names == "1_C12"]
end <- guide_data$ends[guide_data$names == "1_C12"]
chr <- guide_data$seqnames[guide_data$names == "1_C12"]
strand <- guide_data$strands[guide_data$names == "1_C12"]


setwd("C:/Users/yolib/Documents/TFM/Linf_T/Datos/Datos/612-1_C12")
list_bam_full <- list.files(pattern = '*[0-9].bam$')

bam <- "RL384-00022_E24.bam"
name <- "RL384-00022_E24"


ref <- Biostrings::DNAString(reference) 
gd <- GenomicRanges::GRanges(chr, IRanges::IRanges(start + 1, end), strand = strand, genome = "hg38")

names <- str_replace(list_bam_full, ".bam", "")

crispr_set <- readsToTarget(bam, target = gd, reference = ref,
                            target.loc = 30, collapse.pairs = TRUE, names = name)

p <- plotVariants(crispr_set, top.plot = 1,  
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency")) 

