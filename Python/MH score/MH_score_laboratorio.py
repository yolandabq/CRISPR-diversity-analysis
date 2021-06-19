# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 19:49:25 2021

@author: yolib


Calculation of the pattern and MH scores of the experimental data from Lluis Montoliu's lab

"""

#%%
import pandas as pd
#from math import e
import os
import numpy as np 
import glob 
import matplotlib.pyplot as plt
os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\python_scripts')
import MH_functions 
#%%

seqGpr = "aaaacatttccaatgtgaatgcaacagacatttGGCCTGCTACTTTCTGTGTGgggagcgcagtaagttaccttcctttgttcaccctcctccctttgct"

# Firts, we obtain all the possible deletions mediated by MH. 


MH_functions.MH_patterns(seqGpr.upper(), PC = 49, file = '..\\MH\\Gpr143.csv')

# We calculate the MH scores. 


df_Gpr = MH_functions.pattern_score("..\\MH\\Gpr143.csv")

MH_score_Gpr = df_Gpr["Pattern_Score"].sum() 
df_Gpr.to_csv('..\\MH\\Gpr143.csv', sep = "\t")


seqR402Q = "tcctgactctgagtaacccttccctctgtagtaTTTTTGAACAATGGCTGCGAaggcaccgccctcttttggaagtttacccagaagccaatgcacctat"

MH_functions.MH_patterns(seqR402Q.upper(), PC = 49, file = '..\\MH\\R402Q.csv')

df_R402Q = MH_functions.pattern_score("..\\MH\\R402Q.csv")
MH_score_R402Q = df_R402Q["Pattern_Score"].sum()
df_R402Q.to_csv('..\\MH\\R402Q.csv', sep = "\t")


seqS192Y = "tatggatgcattactatgtgtcaagggacacacTGCTTGGGGGCTCTGAAATAtggagggacattgattttgcccatgaagcaccagggtttctgccttg"

MH_functions.MH_patterns(seqS192Y.upper(), PC = 49, file = '..\\MH\\S192Y.csv')

df_S192Y = MH_functions.pattern_score("..\\MH\\S192Y.csv")
MH_score_S192Y = df_S192Y["Pattern_Score"].sum()
df_S192Y.to_csv('..\\MH\\S192Y.csv', sep = "\t")


seqOCA7 = "ttcccccccCTTGCCTTTGTCATTGCAGGTCACTGGAAGGACTGAGTGCATTCAGGAGCCTGGAGGAGCTCATTTTAGACAACAatctgctcggagacga"

MH_functions.MH_patterns(seqOCA7.upper(), PC = 49, file = '..\\MH\\OCA7.csv')

df_OCA7 = MH_functions.pattern_score("..\\MH\\OCA7.csv")
MH_score_OCA7 = df_OCA7["Pattern_Score"].sum()
df_OCA7.to_csv('..\\MH\\OCA7.csv', sep = "\t")


seqH610 = "AATTGAGAAAATGAGTTGTGTTTTTTTTTCTCTTAAAGCACAGGATTTCAGACAGGAGTCTGCTTGTCAAGTGCCTGACGGTGCTGGGAT"

MH_functions.MH_patterns(seqH610.upper(), PC = 49, file = '..\\MH\\H610.csv')

df_H610 = MH_functions.pattern_score("..\\MH\\H610.csv")
MH_score_H610 = df_H610["Pattern_Score"].sum()
df_H610.to_csv('..\\MH\\H610.csv', sep = "\t")



seq476 = "tctcaactcctgattggaaacaatgataacattTGGTGGGTCCCCAATAGCAGTGGcagctcctccaatgtttgtgaagatcacttctgcaatgaggact"

MH_functions.MH_patterns(seq476.upper(), PC = 49, file = '..\\MH\\A476.csv')

df_476 = MH_functions.pattern_score("..\\MH\\A476.csv")
MH_score_476 = df_476["Pattern_Score"].sum()
df_476.to_csv('..\\MH\\A476.csv', sep = "\t")


seqLipo = "acctgccagtgctcaggcaacttcatgggtttcAACTGCGGAAACTCTAAGTTTGGatttgggggcccaaattgtacagagaagcgagtcttgattagaa"

MH_functions.MH_patterns(seqLipo.upper(), PC = 49, file = '..\\MH\\Lipo.csv')

df_Lipo = MH_functions.pattern_score("..\\MH\\Lipo.csv")
MH_score_Lipo = df_Lipo["Pattern_Score"].sum()
df_Lipo.to_csv('..\\MH\\Lipo.csv', sep = "\t")



