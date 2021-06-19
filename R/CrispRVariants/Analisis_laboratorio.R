
# Script to analyze the experiments of the laboratory if Lluis Montoliu with CrispRVariants 

library(CrispRVariants)
library("rtracklayer")
library("GenomicFeatures")
library("gdata")
library(tidyverse)
library("RColorBrewer")

################################################## F0_Gpr143_360+1 ##################################################################

setwd("~/TFM/Alelos_laboratorio/F0_Gpr143/360+1/")

list_bam_full <- list.files(pattern = '*.bam$')
#seq = "aacagacatttGGCCTGCTACTTTCTGTGTGgggagcgcagtaagt"
seq <- "aaaacatttccaatgtgaatgcaacagacatttGGCCTGCTACTTTCTGTGTGgggagcgcagtaagttaccttcctttgttcaccctcctccctttgct"
ref <- Biostrings::DNAString(toupper(seq))
gd <- GenomicRanges::GRanges("chrX", IRanges::IRanges(149219634, 149219733), strand = '+', genome = "mm9")

#list_bam <- c("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9196_S30.bam", "/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam")


#crispr_set <- readsToTarget("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam", target = g, reference = ref,
#target.loc = 20, collapse.pairs = TRUE)

crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 50)


# filtrar el PlotVariants: 


#gdl <- resize(gd, width(gd) + 10, fix = "center")

p <- plotVariants(crispr_set,
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"), 
                  plot.text.size = 2
)

# Mejor hacerlos por separado

plotFreqHeatmap(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)
#plotFreqHeatmap(crispr_set, top.n = 10, type = "proportions", plot.text.size = 2)
plotAlignments(crispr_set, plot.text.size = 2)


################################################## F0_Gpr143 P300L ##################################################################

setwd("~/TFM/Alelos_laboratorio/F0_Gpr143/P300L//")

list_bam_full <- list.files(pattern = '*.bam$')

ref <- Biostrings::DNAString("gactttaatgtattcttttgccttccttcctagggaatactgaatccagcccaaggacttctcttgtctctggccttctatggctggacaggatgcagcc")
gd <- GenomicRanges::GRanges("chrX", IRanges::IRanges(149232888, 149232987), strand = '+', genome = "mm9")

crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 50)

p <- plotVariants(crispr_set, 
                  #plotFreqHeatmap.args = list(top.n = 5, type = "proportions"),
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"), 
                  plot.text.size = 2
)


plotFreqHeatmap(crispr_set, plot.text.size = 2)
#plotFreqHeatmap(crispr_set, top.n = 10, type = "proportions", plot.text.size = 2)
plotAlignments(crispr_set, plot.text.size = 2)



################################################## R402Q Fondo Híbrido ##################################################################

setwd("~/TFM/Alelos_laboratorio/SNAP gene S192Y y R402Q/microinyecc_R402Q_fondo_hibrido/fastq/")

list_bam_full <- list.files(pattern = '*.bam$')

ref <- Biostrings::DNAString("tcctgactctgagtaacccttccctctgtagtaTTTTTGAACAATGGCTGCGAaggcaccgccctcttttggaagtttacccagaagccaatgcacctat")
gd <- GenomicRanges::GRanges("chr7", IRanges::IRanges(94586560, 94586659), strand = '-', genome = "mm9")

#list_bam <- c("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9196_S30.bam", "/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam")


#crispr_set <- readsToTarget("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam", target = g, reference = ref,
#target.loc = 20, collapse.pairs = TRUE)

crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 50)


# filtrar el PlotVariants: 


#gdl <- resize(gd, width(gd) + 10, fix = "center")

p <- plotVariants(crispr_set, #plotAlignments.args = list(top.n = 5), 
                  #plotFreqHeatmap.args = list(top.n = 5, type = "proportions"),
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"), 
                  plot.text.size = 2
)

# Mejor hacerlos por separado

plotFreqHeatmap(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)
#plotFreqHeatmap(crispr_set, top.n = 10, type = "proportions", plot.text.size = 2)
plotAlignments(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)


################################################## R402Q Fondo Puro ##################################################################

setwd("~/TFM/Alelos_laboratorio/SNAP gene S192Y y R402Q/microinyecc_R402Q_fondo_puro/")

list_bam_full <- list.files(pattern = '*.bam$')

ref <- Biostrings::DNAString("tcctgactctgagtaacccttccctctgtagtaTTTTTGAACAATGGCTGCGAaggcaccgccctcttttggaagtttacccagaagccaatgcacctat")
gd <- GenomicRanges::GRanges("chr7", IRanges::IRanges(94586560, 94586659), strand = '-', genome = "mm9")
#list_bam <- c("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9196_S30.bam", "/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam")


#crispr_set <- readsToTarget("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam", target = g, reference = ref,
#target.loc = 20, collapse.pairs = TRUE)

crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 50)


# filtrar el PlotVariants: 


#gdl <- resize(gd, width(gd) + 10, fix = "center")

p <- plotVariants(crispr_set, 
                  #plotFreqHeatmap.args = list(top.n = 5, type = "proportions"),
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"), 
                  plot.text.size = 2
)



plotFreqHeatmap(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)
#plotFreqHeatmap(crispr_set, top.n = 10, type = "proportions", plot.text.size = 2)
plotAlignments(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)



################################################## microinyecc_Doble_R402Q ##################################################################

setwd("~/TFM/Alelos_laboratorio/SNAP gene S192Y y R402Q/microinyeccion_doble/R402Q/")

list_bam_full <- list.files(pattern = '*.bam$')

ref <- Biostrings::DNAString("tcctgactctgagtaacccttccctctgtagtaTTTTTGAACAATGGCTGCGAaggcaccgccctcttttggaagtttacccagaagccaatgcacctat")
gd <- GenomicRanges::GRanges("chr7", IRanges::IRanges(94586560, 94586659), strand = '-', genome = "mm9")
#list_bam <- c("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9196_S30.bam", "/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam")


#crispr_set <- readsToTarget("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam", target = g, reference = ref,
#target.loc = 20, collapse.pairs = TRUE)

crispr_set <- readsToTarget(list_bam_full[c(3:14, 16, 17, 18, 19, 20, 21)], target = gd, reference = ref,
                            target.loc = 66, chimeras = "merge")


# filtrar el PlotVariants: 


#gdl <- resize(gd, width(gd) + 10, fix = "center")

p <- plotVariants(crispr_set, 
                  #plotFreqHeatmap.args = list(top.n = 5, type = "proportions"),
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"), 
                  plot.text.size = 2
)

ch <- getChimeras(crispr_set, sample =  "B9079.bam")
plotChimeras(ch)



plotFreqHeatmap(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)
#plotFreqHeatmap(crispr_set, top.n = 10, type = "proportions", plot.text.size = 2)
plotAlignments(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)


################################################## microinyecc_Doble_S192Y ##################################################################

setwd("~/TFM/Alelos_laboratorio/SNAP gene S192Y y R402Q/microinyeccion_doble/S192Y/")

list_bam_full <- list.files(pattern = '*.bam$')

ref <- Biostrings::DNAString("tatggatgcattactatgtgtcaagggacacacTGCTTGGGGGCTCTGAAATAtggagggacattgattttgcccatgaagcaccagggtttctgccttg")
gd <- GenomicRanges::GRanges("chr7", IRanges::IRanges(94641232, 94641328), strand = '-', genome = "mm9")
#list_bam <- c("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9196_S30.bam", "/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam")


#crispr_set <- readsToTarget("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam", target = g, reference = ref,
#target.loc = 20, collapse.pairs = TRUE)

crispr_set <- readsToTarget(list_bam_full[c(1:10, 12:31)], target = gd, reference = ref,
                            target.loc = 50)

# 1:10, 12:20, 20:25, 
# filtrar el PlotVariants: 


#gdl <- resize(gd, width(gd) + 10, fix = "center")

p <- plotVariants(crispr_set, 
                  #plotFreqHeatmap.args = list(top.n = 5, type = "proportions"),
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"), 
                  plot.text.size = 2
)

ch <- getChimeras(crispr_set, sample =  "B9079.bam")
plotChimeras(ch)



plotFreqHeatmap(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)
#plotFreqHeatmap(crispr_set, top.n = 10, type = "proportions", plot.text.size = 2)
plotAlignments(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)


################################################## OCA7 ##################################################################

setwd("~/TFM/Alelos_laboratorio/OCA7/lith/Fastq/fastq_merged/")

list_bam_full <- list.files(pattern = '*.bam$')

ref <- Biostrings::DNAString("ttcccccccCTTGCCTTTGTCATTGCAGGTCACTGGAAGGACTGAGTGCATTCAGGAGCCTGGAGGAGCTCATTTTAGACAACAatctgctcggagacga")
gd <- GenomicRanges::GRanges("chr14", IRanges::IRanges(23397042, 23397141), strand = '+', genome = "mm9")
#list_bam <- c("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9196_S30.bam", "/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam")


#crispr_set <- readsToTarget("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam", target = g, reference = ref,
#target.loc = 20, collapse.pairs = TRUE)

crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 50)

# 1:10, 12:20, 20:25, 
# filtrar el PlotVariants: 


#gdl <- resize(gd, width(gd) + 10, fix = "center")

p <- plotVariants(crispr_set, 
                  #plotFreqHeatmap.args = list(top.n = 5, type = "proportions"),
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"), 
                  plot.text.size = 2
)


##################################################  OCA2_H610 ##################################################################

# OCA2_H610

setwd("~/TFM/Alelos_laboratorio/OCA2_H610")

list_bam_full <- list.files(pattern = '*.bam$')

#names <- c(30:65)
a <- 'AATT GAGAAAATGA GTTGTGTTTT TTTTTCTCTT AAAGCaCAGg  63587334ATTTCAGACA GGAGTCTGCT TGTCAAGTGC CTGACGGTGC tgggat'

a2 <- gsub("[0-9]","",a)
a3 <- gsub(" ", "", a2)
a4 <- gsub("\n","", a3)

ref <- Biostrings::DNAString(a4)

gd <- GenomicRanges::GRanges("chr7", IRanges::IRanges(63587291, 63587380), strand = '+', genome = "mm9")


crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 50, collapse.pairs = TRUE, names = list_bam_full)

p <- plotVariants(crispr_set, 
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency"),
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"), 
                  plot.text.size = 2
)


plotFreqHeatmap(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)
#plotFreqHeatmap(crispr_set, top.n = 10, type = "proportions", plot.text.size = 2)
plotAlignments(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)


##################################################  OCA2_476 ##################################################################


setwd("~/TFM/Alelos_laboratorio/OCA2 476")

list_bam_full <- list.files(pattern = '*.bam$')

names <- c(30:65)

ref <- Biostrings::DNAString("tctcaactcctgattggaaacaatgataacattTGGTGGGTCCCCAATAGCAGTGGcagctcctccaatgtttgtgaagatcacttctgcaatgaggact")
gd <- GenomicRanges::GRanges("chr7", IRanges::IRanges(63580057, 63580156), strand = '-', genome = "mm9")

#list_bam <- c("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9196_S30.bam", "/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam")


#crispr_set <- readsToTarget("/home/ubuntu/Desktop/TFM/Alelos_laboratorio/OCA2_476/B9197_S31.bam", target = g, reference = ref,
#target.loc = 20, collapse.pairs = TRUE)

crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 50, collapse.pairs = TRUE, names = names)


# filtrar el PlotVariants: 


#gdl <- resize(gd, width(gd) + 10, fix = "center")

p <- plotVariants(crispr_set, plotAlignments.args = list(top.n = 5), 
                  plotFreqHeatmap.args = list(top.n = 5, type = "counts"),
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"), 
                  plot.text.size = 2
)

# Mejor hacerlos por separado

plotFreqHeatmap(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2,  header = "efficiency")
#plotFreqHeatmap(crispr_set, top.n = 10, type = "proportions", plot.text.size = 2)
plotAlignments(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)


##################################################  Tyr (lipos) ##################################################################

# Now I am going to analyze the data from the lipos experiments. N2a cells are mutated
# and we want to restore the normal sequence.
# The unmuted sequence, the one normally found in mouse cells, is: "AACTGCGGAAACTGTAAATTTGG"
# The one we use as a guide in this case, which will be attached to the mutated sequence, is "AACTGCGGAAACTCTAAGTTTGG"
# If we look for it in BLAT, it cannot find it, because it has the mutation.


setwd("~/TFM/Alelos_laboratorio/Experimentos_Almudena/bam/Lipoparticulas/")

list_bam_full <- list.files(pattern = '*.bam$')
names <- str_replace(list_bam_full, ".bam", "")
#names <- c(30:65)


ref <- Biostrings::DNAString("acctgccagtgctcaggcaacttcatgggtttcAACTGCGGAAACTCTAAGTTTGGatttgggggcccaaattgtacagagaagcgagtcttgattagaa") # voy a utilizar la mutada
gd <- GenomicRanges::GRanges("chr7", IRanges::IRanges(94641500, 94641599), strand = '-', genome = "mm9")


crispr_set <- readsToTarget(list_bam_full, target = gd, reference = ref,
                            target.loc = 50, collapse.pairs = TRUE, names = names)


common_alleles <- as.data.frame(variantCounts(crispr_set, result = "proportions"),
                                row.names = rownames(variantCounts(crispr_set, result = "proportions")))

common_alleles['Alleles'] = rownames(common_alleles) # para guardar el nombre de la columna 


file_path <- str_c("~/TFM/Alelos_laboratorio/Experimentos_Almudena/bam/Lipoparticulas/tyr103_rna.csv") # para guardar el csv
#write.csv(common_alleles, file_path, row.names = FALSE)


#crispr_set_opposite <- readsToTarget(list_bam_full, target = gd, reference = ref,
#target.loc = 17, collapse.pairs = TRUE, orientation = "opposite")

grp <- c(1,2,1,2)
col_grp <- c('red', 'blue')

p <- plotVariants(crispr_set, top.plot = 1, plot.title =  "Alelos resultantes segÃºn la lipoproteÃ­na utilizada", 
                  plotFreqHeatmap.args = list( type = "proportions", header = "efficiency",
                                               group = grp, group.colours = col_grp)) # ahora el header es la eficiencia 

plotFreqHeatmap(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)
#plotFreqHeatmap(crispr_set, top.n = 10, type = "proportions", plot.text.size = 2)
plotAlignments(crispr_set, min.freq = 10, type = "proportions", plot.text.size = 2)


crispr_set$mutationEfficiency()

prop_crispr_set <- prop.table(crispr_set$cigar_freqs, 2)*100

df_prop_crispr_set <- as.data.frame(prop_crispr_set, row.names = rownames(prop_crispr_set))

# We prepare the data for the pie charts. 

#  L3000_DNA_S21

library(dplyr)
library("devtools")

L3000_DNA_S33 <- df_prop_crispr_set %>% select(L3000_DNA_S33)
L3000_DNA_S33$alleles <- rownames(L3000_DNA_S33)
L3000_DNA_S33$alleles[L3000_DNA_S33$alleles == 'no variant'] <- 'Wild Type'

df_L3000_DNA_S33_types <- L3000_DNA_S33[L3000_DNA_S33$L3000_DNA_S33>1,]
df_L3000_DNA_S33_types <- rbind(df_L3000_DNA_S33_types, c(100-sum(df_L3000_DNA_S33_types$L3000_DNA_S33), "Others"))
colnames(df_L3000_DNA_S33_types) <- c("porcentage", "alleles")
n_alleles <- length(df_L3000_DNA_S33_types$porcentage)

# L3000

#prop_crispr_set <- prop.table(crispr_set$cigar_freqs, 2)*100

#df_prop_crispr_set <- as.data.frame(prop_crispr_set, row.names = rownames(prop_crispr_set))

df_L300 <- df_prop_crispr_set %>% select(L3000_S21)
df_L300$alleles <- rownames(df_L300)
df_L300$alleles[df_L300$alleles == 'no variant'] <- 'Wild Type'

df_L300_types <- df_L300[df_L300$L3000_S21>1,]
df_L300_types <- rbind(df_L300_types, c(100-sum(df_L300_types$L3000_S21), "Others"))
colnames(df_L300_types) <- c("porcentage", "alleles")
n_alleles <- length(df_L300_types$porcentage)


# LMAX_S22.bam


df_LMAX_S22 <- df_prop_crispr_set %>% select(LMAX_S22)
df_LMAX_S22$alleles <- rownames(df_LMAX_S22)
df_LMAX_S22$alleles[df_LMAX_S22$alleles == 'no variant'] <- 'Wild Type'

df_LMAX_S22_types <- df_LMAX_S22[df_LMAX_S22$LMAX_S22>1,]
df_LMAX_S22_types <- rbind(df_LMAX_S22_types, c(100-sum(df_LMAX_S22_types$LMAX_S22), "Others"))
colnames(df_LMAX_S22_types) <- c("porcentage", "alleles")
n_alleles <- length(df_LMAX_S22_types$porcentage)


# LMAX_DNA_S34


df_LMAX_DNA_S34 <- df_prop_crispr_set %>% select(LMAX_DNA_S34)
#df_LMAX_DNA_S34 <- sort(df_LMAX_DNA_S34)
df_LMAX_DNA_S34$alleles <- rownames(df_LMAX_DNA_S34)
df_LMAX_DNA_S34$alleles[df_LMAX_DNA_S34$alleles == 'no variant'] <- 'Wild Type'

df_LMAX_DNA_S34_types <- df_LMAX_DNA_S34[df_LMAX_DNA_S34$LMAX_DNA_S34>1,]
df_LMAX_DNA_S34_types <- rbind(df_LMAX_DNA_S34_types, c(100-sum(df_LMAX_DNA_S34_types$LMAX_DNA_S34), "Others"))
colnames(df_LMAX_DNA_S34_types) <- c("porcentage", "alleles")
n_alleles <- length(df_LMAX_DNA_S34_types$porcentage)

# We prepare the colors ----------

alleles_list <- list()
alleles_list <- c(df_L3000_DNA_S33_types$alleles, df_LMAX_DNA_S34_types$alleles, df_LMAX_S22_types$alleles,
                  df_L300_types$alleles)

alleles_unique <- unique(alleles_list)

col_alleles <- data.frame(
  alleles = alleles_unique,  
  col = brewer.pal(n = length(alleles_unique), name = 'Set3'),
  stringsAsFactors = FALSE)


#We do the plots ----------------

df_L3000_DNA_S33_types$col <- col_alleles$col[match(df_L3000_DNA_S33_types$alleles, col_alleles$alleles)]
#df_L300_DNA_S21_types$col <- as.factor(df_L300_DNA_S21_types$col)
df_L3000_DNA_S33_types <- df_L3000_DNA_S33_types[order(as.numeric(df_L3000_DNA_S33_types$porcentage), decreasing = TRUE), ]


ggplot(df_L3000_DNA_S33_types,aes(x="",y=as.numeric(porcentage), fill=alleles))+
  geom_bar(stat = "identity",color="white")+
  coord_polar(theta="y") + 
  scale_fill_manual(breaks = df_L3000_DNA_S33_types$alleles, values = df_L3000_DNA_S33_types$col)+
  #values=col_alleles$col[col_alleles$alleles == df_L300_DNA_S21_types$alleles])+
  #scale_fill_manual(values=brewer.pal(n = n_alleles, name = 'Paired')) + 
  geom_text(aes(label=paste(round(as.numeric(porcentage), 2), "%"), x = 1.3), 
            position=position_stack(vjust=0.5), color="black",size=3) +
  labs(y = NULL,  x = NULL) + 
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  labs(y = NULL) + 
  ggtitle("Distribución de alelos en L3000_DNA_S1")


#display.brewer.pal(n = n_alleles, name = 'RdBu')

df_L300_types$col <- col_alleles$col[match(df_L300_types$alleles, col_alleles$alleles)]
#df_L300_types$alleles <- as.factor(df_L300_types$alleles)
df_L300_types <- df_L300_types[order(as.numeric(df_L300_types$porcentage), decreasing = TRUE), ]


ggplot(df_L300_types,aes(x="",y=as.numeric(porcentage), fill=alleles))+
  geom_bar(stat = "identity",color="white")+
  coord_polar(theta="y", start = 0) + 
  scale_fill_manual(breaks = df_L300_types$alleles,values=df_L300_types$col)+
  #scale_fill_manual(values=brewer.pal(n = n_alleles, name = 'Paired')) + 
  geom_text(aes(label=paste(round(as.numeric(porcentage), 2), "%"), x = 1.3), 
            position=position_stack(vjust=0.5), color="black",size=3) +
  labs(y = NULL,  x = NULL) + 
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + labs(y = NULL) + 
  ggtitle("Distribución de alelos en L3000")


#display.brewer.pal(n = n_alleles, name = 'RdBu')


df_LMAX_S22_types$col <- col_alleles$col[match(df_LMAX_S22_types$alleles, col_alleles$alleles)]
df_LMAX_S22_types <- df_LMAX_S22_types[order(as.numeric(df_LMAX_S22_types$porcentage), decreasing = TRUE), ]


ggplot(df_LMAX_S22_types,aes(x="",y=as.numeric(porcentage), fill=alleles))+
  geom_bar(stat = "identity",color="white")+
  coord_polar(theta="y") + 
  scale_fill_manual(breaks = df_LMAX_S22_types$alleles,values=df_LMAX_S22_types$col)+
  #scale_fill_manual(values=brewer.pal(n = n_alleles, name = 'Paired')) + 
  geom_text(aes(label=paste(round(as.numeric(porcentage), 2), "%"), x = 1.3), 
            position=position_stack(vjust=0.5), color="black",size=3) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  labs(y = NULL,  x = NULL) + 
  ggtitle("Distribución de alelos en LMAX_S22")


#df_LMAX_DNA_S34_types <- df_LMAX_DNA_S34_types[order(as.numeric(df_LMAX_DNA_S34_types$porcentage), decreasing = TRUE), ]
#display.brewer.pal(n = n_alleles, name = 'RdBu')

df_LMAX_DNA_S34_types$col <- col_alleles$col[match(df_LMAX_DNA_S34_types$alleles, col_alleles$alleles)]
df_LMAX_DNA_S34_types <- df_LMAX_DNA_S34_types[order(as.numeric(df_LMAX_DNA_S34_types$porcentage), decreasing = TRUE), ]


ggplot(df_LMAX_DNA_S34_types,aes(x="",y=as.numeric(porcentage), fill=alleles))+
  geom_bar(stat = "identity",color="white")+
  coord_polar(theta="y") + 
  scale_fill_manual(breaks = df_LMAX_DNA_S34_types$alleles,values=df_LMAX_DNA_S34_types$col)+
  geom_text(aes(label=paste(round(as.numeric(porcentage), 2), "%"), x = 1.3), 
            position=position_stack(vjust=0.5), color="black",size=3) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  labs(y = NULL, x = NULL) + 
  ggtitle("Distribución de alelos en LMAX-DNA_S34")





