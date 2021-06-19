# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 11:14:44 2021

@author: yolib

Analysis of deletions larger than 3 bp mediated by MH. 

"""

#%%

import os
import pandas as pd
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
os.chdir('C:\\Users\\yolib\\OneDrive\\Documentos\\TFM\\Documento TFM\\Codigo\\Python\\Funciones')
from MH_functions import MH_detection


#%% Verify if deletions larger than 3 bp are mediated by MH

""" 
Analysis of deletions that overlap with the cleavage site. 
If some mutation doesn't overlap with the DSB, I try to modify the alignment. 

"""
"""
os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\')

guides_seq = pd.read_csv("guides_sequence_2.csv")

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results')
"""

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\')

guides_seq = pd.read_csv("guides_sequence_2_no_replic.csv")

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results_no_replic')

guides_files = os.listdir()

MH_mut_list = []
NO_MH_mut_list = []
MH_seq_list = []

#data = pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results\\' + file)

for guide in guides_files: 
    
    data = pd.read_csv(guide)
    data_sorted = data.sort_values(list(data.columns)[:-1], ascending = False)
    #data_0_5_freq_index = data_sorted[data_sorted.loc[:, list(data_sorted.columns)[:-1]] >= 0.5 ].dropna(how='all').index
    #data_sorted = data_sorted.loc[data_0_5_freq_index, :]
    data_sorted = data_sorted[data_sorted["Alleles"] != "no variant"]
    data_sorted = data_sorted.iloc[:25, :]
    
    seq = guides_seq[guides_seq["Guide_ID"] == int(guide.split("_")[0])]["Sequence"].iloc[0]
    
    alleles = [cigar for cigar in data_sorted["Alleles"] if cigar[-1]=="D" and "," not in cigar]
    alleles_mayor_3 = [cigar for cigar in alleles if int(cigar.split(":")[1][:-1]) >=3 ]
        # We set < 60 to avoid sequences with long deletions that go out of the 60 bp that we have
    for allele in alleles_mayor_3: 
        
        in_PC = True
        MH_sequence = ''
        
        if int(allele.split(":")[0]) < 0 and (int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) >=0 and (30 + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= 60: 
            MH_sequence = MH_detection(seq, allele)
            
        elif int(allele.split(":")[0]) == 1 and (30 + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= 60: # if the starting position is +1, the deletion overlap with the cleavage site
            
            MH_sequence = MH_detection(seq, allele)
            
        elif int(allele.split(":")[0]) > 1 and (30 + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= 60: 
        # If the starting position of the indel is after the cleavage site, we try to realign the sequence and we verify again if now the deletion overlap with the DSB. 
            pos_inic_indel = 30 + int(allele.split(":")[0]) - 2
            alt_alin = True
            alt_pos = 0
            
            while alt_alin:
                
                if seq[pos_inic_indel] == seq[pos_inic_indel + int(allele.split(":")[1][:-1])]: 
                    alt_pos = alt_pos + 1
                    pos_inic_indel = pos_inic_indel - 1
        
                else: 
                    alt_alin = False
                    
            alt_allele_cigar = int(allele.split(":")[0]) - alt_pos
            alt_allele_cigar_str = str(alt_allele_cigar) + ":" + allele.split(":")[1][:-1] + "D"
            
            if int(alt_allele_cigar_str.split(":")[0]) < 0 and (int(alt_allele_cigar_str.split(":")[0]) + int(alt_allele_cigar_str.split(":")[1][:-1])) >=0: 
                MH_sequence = MH_detection(seq, alt_allele_cigar_str)
                
            elif int(alt_allele_cigar_str.split(":")[0]) == 1: 
                
                MH_sequence = MH_detection(seq, alt_allele_cigar_str)
            else: 
                
                in_PC = False
                
        elif int(allele.split(":")[0]) < 0 and (int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) < 0 and (30 + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= 60 and int(allele.split(":")[0]) >= -60: 
              # If the deletion start before the cleavage site but it also ends before the cleavage site, we try to realign the sequence. 
            pos_inic_indel = 30 + int(allele.split(":")[0]) - 1
            alt_alin = True
            alt_pos = 0
            while alt_alin:
                if seq[pos_inic_indel + 1] == seq[pos_inic_indel + 1 + int(allele.split(":")[1][:-1])]: 
                    alt_pos = alt_pos + 1
                    pos_inic_indel = pos_inic_indel + 1
        
                else: 
                    alt_alin = False
                    
            alt_allele_cigar = int(allele.split(":")[0]) + alt_pos
            alt_allele_cigar_str = str(alt_allele_cigar) + ":" + allele.split(":")[1][:-1] + "D"
            
            if int(alt_allele_cigar_str.split(":")[0]) < 0 and (int(alt_allele_cigar_str.split(":")[0]) + int(alt_allele_cigar_str.split(":")[1][:-1])) >=0: 
                MH_sequence = MH_detection(seq, alt_allele_cigar_str)
                
            elif int(alt_allele_cigar_str.split(":")[0]) == 1: 
                
                MH_sequence = MH_detection(seq, alt_allele_cigar_str)
    
            else: 
                
                in_PC = False
                
        if MH_sequence == '' and in_PC == True: 
            
            NO_MH_mut_list.append(allele)
            
        elif MH_sequence != '' and in_PC == True: 
            MH_mut_list.append(allele)
            MH_seq_list.append(MH_sequence)
            
            

#%%

print("% of deletions >= 3 bp with MH (top 25 alleles): " + str(round(len(MH_mut_list)/(len(MH_mut_list) + len(NO_MH_mut_list)), 2)))
print("% of deletions >= 3 bp without MH (top 25 alleles): " + str(round(len(NO_MH_mut_list)/(len(MH_mut_list) + len(NO_MH_mut_list)), 2)))



#%%

colors = sns.color_palette("Set2", 2)


plt.pie([len(MH_mut_list), len(NO_MH_mut_list)], colors = colors, 
        labels=["Deletions > 3  with MH","Deletions > 3 without MH"], autopct='%1.1f%%',
        startangle=90)

#plt.title("Distribución de MH en deleciones > 3 \n (alelos con frecuencia> 0.5)")
plt.title("MH events in deletions > 3 bp \n (top 25 alleles)")

#%% Size of the MH

MH_size = [len(x) for x in MH_seq_list]
#del_size = [int(x.split(":")[1][:-1]) for x in MH_mut_list]

from collections import Counter

a = dict(Counter(MH_size))
print(a)

plt.figure(figsize=(9, 6))

sns.histplot(MH_size ,discrete=True, stat = "probability")
#plt.title("Distribución del tamaño de MH en deleciones > 2 \n (alelos con frecuencia> 0.5)")
plt.title("Size of MH associated with deletions events \n (top 25 alleles)", fontsize = 15)
plt.xlabel("Size of the MH", size = 12)

MH_size = np.array(MH_size)
print("% of MH with length 2: " + str(100*len(MH_size[MH_size == 2]) / len(MH_size)))
print("% of MH with length 3: " + str(100*len(MH_size[MH_size == 3]) / len(MH_size)))
print("% of MH with length 4: " + str(100*len(MH_size[MH_size == 4]) / len(MH_size)))
print("% of MH with length 5: " + str(100*len(MH_size[MH_size == 5]) / len(MH_size)))

print("% of MH with length > 5: " + str(100*len(MH_size[MH_size > 5]) / len(MH_size)))


#%% All the alleles (not only top 25)

"""
os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\')

guides_seq = pd.read_csv("guides_sequence_2.csv")

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results')
"""

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\')

guides_seq = pd.read_csv("guides_sequence_2_no_replic.csv")

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results_no_replic')

guides_files = os.listdir()

MH_mut_list = []
NO_MH_mut_list = []
MH_seq_list = []

#data = pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results\\' + file)

for guide in guides_files: 
    
    data = pd.read_csv(guide)
    data_sorted = data.sort_values(list(data.columns)[:-1], ascending = False)
    data_0_5_freq_index = data_sorted[data_sorted.loc[:, list(data_sorted.columns)[:-1]] >= 0 ].dropna(how='all').index
    data_sorted = data_sorted.loc[data_0_5_freq_index, :]
    seq = guides_seq[guides_seq["Guide_ID"] == int(guide.split("_")[0])]["Sequence"].iloc[0]
    
    alleles = [cigar for cigar in data_sorted["Alleles"] if cigar[-1]=="D" and "," not in cigar]
    alleles_mayor_3 = [cigar for cigar in alleles if int(cigar.split(":")[1][:-1]) >=3 ]
    # lo de menor de 60 lo pongo para evitar que haya errores de secuencias con deleciones muy grandes que se salgan de la secuencia de 60 pb que tenemos
    for allele in alleles_mayor_3: 
        
        in_PC = True
        MH_sequence = ''
        
        if int(allele.split(":")[0]) < 0 and (int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) >=0 and (30 + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= 60: 
            MH_sequence = MH_detection(seq, allele)
            
        elif int(allele.split(":")[0]) == 1 and (30 + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= 60: 
            
            MH_sequence = MH_detection(seq, allele)
            
        elif int(allele.split(":")[0]) > 1 and (30 + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= 60:  
       
            pos_inic_indel = 30 + int(allele.split(":")[0]) - 2
            alt_alin = True
            alt_pos = 0
            
            while alt_alin:
                
                if seq[pos_inic_indel] == seq[pos_inic_indel + int(allele.split(":")[1][:-1])]: 
                    alt_pos = alt_pos + 1
                    pos_inic_indel = pos_inic_indel - 1
        
                else: 
                    alt_alin = False
                    
            alt_allele_cigar = int(allele.split(":")[0]) - alt_pos
            alt_allele_cigar_str = str(alt_allele_cigar) + ":" + allele.split(":")[1][:-1] + "D"
            
            if int(alt_allele_cigar_str.split(":")[0]) < 0 and (int(alt_allele_cigar_str.split(":")[0]) + int(alt_allele_cigar_str.split(":")[1][:-1])) >=0: 
                MH_sequence = MH_detection(seq, alt_allele_cigar_str)
                
            elif int(alt_allele_cigar_str.split(":")[0]) == 1: 
                
                MH_sequence = MH_detection(seq, alt_allele_cigar_str)
            else: 
                
                in_PC = False
                
        elif int(allele.split(":")[0]) < 0 and (int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) < 0 and (30 + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= 60 and int(allele.split(":")[0]) >= -60: 
         
            pos_inic_indel = 30 + int(allele.split(":")[0]) - 1
            alt_alin = True
            alt_pos = 0
            while alt_alin:
                if seq[pos_inic_indel + 1] == seq[pos_inic_indel + 1 + int(allele.split(":")[1][:-1])]: 
                    alt_pos = alt_pos + 1
                    pos_inic_indel = pos_inic_indel + 1
        
                else: 
                    alt_alin = False
                    
            alt_allele_cigar = int(allele.split(":")[0]) + alt_pos
            alt_allele_cigar_str = str(alt_allele_cigar) + ":" + allele.split(":")[1][:-1] + "D"
            
            if int(alt_allele_cigar_str.split(":")[0]) < 0 and (int(alt_allele_cigar_str.split(":")[0]) + int(alt_allele_cigar_str.split(":")[1][:-1])) >=0: 
                MH_sequence = MH_detection(seq, alt_allele_cigar_str)
                
            elif int(alt_allele_cigar_str.split(":")[0]) == 1: 
                
                MH_sequence = MH_detection(seq, alt_allele_cigar_str)
    
            else: 
                
                in_PC = False
                
        if MH_sequence == '' and in_PC == True: 
            
            NO_MH_mut_list.append(allele)
            
        elif MH_sequence != '' and in_PC == True: 
            MH_mut_list.append(allele)
            MH_seq_list.append(MH_sequence)
            

#%%


print("% of deletions >= 3 bp with MH: " + str(round(len(MH_mut_list)/(len(MH_mut_list) + len(NO_MH_mut_list)), 2)))
print("% of deletions >= 3 bp without MH: " + str(round(len(NO_MH_mut_list)/(len(MH_mut_list) + len(NO_MH_mut_list)), 2)))




#%% Size of the MH


MH_size = [len(x) for x in MH_seq_list]


from collections import Counter

a = dict(Counter(MH_size))
print(a)

sns.histplot(MH_size ,discrete=True)
plt.title("Distribución del tamaño de MH en deleciones > 2 \n (todos los alelos)")
plt.xlabel("Tamaño de la MH")
plt.show()

#%%

colors = sns.color_palette("Set2", 2)


plt.pie([len(MH_mut_list), len(NO_MH_mut_list)], colors = colors, 
        labels=["Deleciones > 2 con MH","Deleciones > 2 sin MH"], autopct='%1.1f%%',
        startangle=90)

plt.title("Distribución de MH en deleciones > 2 \n (todos los alelos)")
plt.show()