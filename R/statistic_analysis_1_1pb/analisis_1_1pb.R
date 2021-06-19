# Statistical analysis of insertions 1: 1I.

# T test to see if the nucleotide inserted at position 1 is related to the nucleotide at position -1.

library(readr)
ins_nt_menos1 <- read_csv("C:/Users/yolib/Documents/TFM/Linf_T/Datos/analisis_ins/ins_nt_menos1.csv")

ins_nt_menos1$nt_menos_1 <- as.factor(ins_nt_menos1$nt_menos_1)
ins_nt_menos1$nt_inserted <- as.factor(ins_nt_menos1$nt_inserted)


library(dbplyr)
library(dplyr)



library("ggpubr")
library("rstatix")

ins_nt_menos1 <- read_csv("C:/Users/yolib/Documents/TFM/Linf_T/Datos/analisis_ins/ins_nt_menos1.csv")

ins_nt_menos1$nt_menos_1 <- as.factor(ins_nt_menos1$nt_menos_1)
ins_nt_menos1$nt_inserted <- as.factor(ins_nt_menos1$nt_inserted)

a <- ins_nt_menos1[ins_nt_menos1["nt_menos_1"] == "A",]
a2 <- a[a["nt_inserted"] == "A", ]

qqnorm(a$frequency, pch = 1, frame = FALSE)
qqline(a$frequency, col = "steelblue", lwd = 2)

stat.test <- ins_nt_menos1 %>%
  group_by(nt_menos_1) %>%
  #t_test(frequency ~ nt_inserted) %>%
  wilcox_test(frequency ~ nt_inserted) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test

stat.test <- stat.test %>% add_xy_position(x = "nt_inserted")
bxp <- ggboxplot(ins_nt_menos1, x = "nt_inserted", y = "frequency", fill = "#00AFBB",
                 facet.by = "nt_menos_1")
bxp + labs( x = "Nucleotide Inserted") + 

  stat_pvalue_manual(stat.test, label = "p.adj.signif") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))



