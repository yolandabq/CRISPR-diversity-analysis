# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 17:40:30 2021

@author: yolib

Script to get a list of sequences and their "weights", to use this fila as input in kpLogo. 

"""



#%%

import os
import pandas as pd
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
import glob
from collections import Counter
import math
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
os.chdir('C:\\Users\\yolib\\OneDrive\\Documentos\\TFM\\Documento TFM\\Codigo\\Python\\Funciones')
import utils_functions
import entropy_functions as enp


#%%


os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos')
guides_data_no_replic = pd.read_csv("guides_sequence_2_no_replic.csv")
guides_replic_data = pd.read_csv("guides_sequence_2.csv")

frames = [guides_replic_data, guides_data_no_replic]
guides_data = pd.concat(frames)


os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results_no_replic_ins')
guides_ins_files = os.listdir()

#ins_data = pd.read_csv(guides_ins_files[216])


os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results_no_replic')

guides_files = os.listdir()

guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))


#%% We compute the cutting efficiency

seq_dict = {}
seq_list = {}
seq_list["A"] = []
seq_list["C"] = []
seq_list["T"] = []
seq_list["G"] = []

with open('..\\Logos\\file_sequence_efficiency.txt', 'w') as f:
    
    for g_id in guides_id:
        #print(g_id)
        #print(glob.glob(str(g_id) + "_*"))
        files = glob.glob(str(g_id) + "_*")
        #sequence = list(guides_data [guides_data["Guide_ID"] == g_id]["Sequence"])[0]
        porc_corte = []
        for file in files: 
            
            df_results = pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results_no_replic\\' + file)
            porc_corte.append(100-np.mean(list(df_results[df_results["Alleles"] == "no variant"].iloc[0, :-1])))
            
        corte_medio = (np.mean(porc_corte)/100)
        #sequence = (list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0])[24:34].upper()
        #[24:34] si solo quiero la secuencia cerca del punto de corte. 
        
        sequence = (list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0]).upper()

        
        f.write(sequence + "\t" + str(corte_medio) + "\n")
        #seq_dict[sequence] = corte_medio
        seq_list[sequence[8].upper()].append(corte_medio)

#%%
import operator
seq_dict_sort = sorted(seq_dict.items(), key=operator.itemgetter(1), reverse=True)

with open('..\\Logos\\file_sequence_efficiency_ranked.txt', 'w') as f:
    
    [f.write(seq[0] + '\n') for seq in seq_dict_sort]


#%% % of insertions

#guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))

seq_list = {} # We save the sequences depending on the X position you choose
seq_list["A"] = []
seq_list["C"] = []
seq_list["T"] = []
seq_list["G"] = []

guides_nt_1 = {} # We save the sequences depending on the nucleotide in -1
guides_nt_1["A"] = []
guides_nt_1["C"] = []
guides_nt_1["T"] = []
guides_nt_1["G"] = []

guides_nt_3 = {} # We save the sequences depending on the nucleotide in +3
guides_nt_3["A"] = []
guides_nt_3["C"] = []
guides_nt_3["T"] = []
guides_nt_3["G"] = []


df_guias_nt = pd.DataFrame(index = guides_id, columns=["-2", "-1", "+1", "+3"])

guide_insertion = {}

with open('..\\file_sequence_insertion.txt', 'w') as f:
    
    for carpeta in ["results", "results_no_replic"]:
        
        os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta)
        guides_files = os.listdir()
        guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))

        for g_id in guides_id:
            #print(g_id)
            #print(glob.glob(str(g_id) + "_*"))
            files = glob.glob(str(g_id) + "_*")
            #sequence = list(guides_data [guides_data["Guide_ID"] == g_id]["Sequence"])[0]
            #porc_corte_list = []
            porc_ins_list = []
            
            for file in files: 
                
                df_results = pd.read_csv(file)
                #df_results = df_results.sort_values(list(df_results.columns)[0:-1],ascending=False)
                
                if "no variant" in list(df_results["Alleles"]):
                    porc_corte = (100-np.mean(list(df_results[df_results["Alleles"] == "no variant"].iloc[0, :-1])))
                    #porc_corte_list.append(porc_corte)
                
                else: 
                    porc_corte = 100
                    
                porc_ins = 0
                norm = 1/porc_corte
                
                for row in list(df_results.index):
                    #print(row)
                    
                    if ',' not in df_results.loc[row, list(df_results.columns)[-1]] and df_results.loc[row, list(df_results.columns)[-1]][-1] == "I":
                        porc_ins = porc_ins + norm*np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]]))
                        
                porc_ins_list.append(porc_ins)
                #print("ins " + str(porc_in))
         
                #print("****")
                
            
            porc_ins_medio = np.mean(porc_ins_list)
            #corte_medio = (np.mean(porc_corte)/100)
            sequence = (list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0]).upper()
            #sequence = (list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0])[24:34].upper()
            #seq_list[sequence[8].upper()].append(porc_ins_medio)
            seq_list[sequence[8].upper()].append(porc_ins_medio) # La pos 5 es la "-1" y la 8 la "+3". Guardo el % de inserción
            
            df_guias_nt["Guide_ID"] = df_guias_nt.index
            df_guias_nt.loc[g_id, :-1] = [sequence[28].upper(),sequence[29].upper(),sequence[30].upper(),sequence[32].upper()]
            
            
            guides_nt_1[sequence[29].upper()].append(g_id) # guias que tienen "A","T","G" o "C" en el nt -1
            guides_nt_3[sequence[32].upper()].append(g_id)
            
            guide_insertion[g_id] = porc_ins_medio
            f.write(sequence + "\t" + str(porc_ins_medio) + "\n")
            
        df_guias_nt = df_guias_nt.sort_values(by = "Guide_ID")


#%%
# We check for insert differences with different combinations of nt in positions -1 and +3
# For example if there is T at -1 and if there is A at +3 it seems to increase the % of insertions, so we see what
# occurs when both things happen at the same time, and separately.
# Position -2 also seems important (C increases insertion and T decreases it)

#inter = set(guides_nt_1["T"]).intersection(set(guides_nt_3["A"]))
#only_T_1 = set(guides_nt_1["T"]).difference(set(guides_nt_3["A"]))

#only_A_3 = set(guides_nt_3["A"]).difference(set(guides_nt_1["T"]))

inter = set(list(df_guias_nt[df_guias_nt["-1"] == "T"].index)).intersection(set(list(df_guias_nt[df_guias_nt["+3"] == "A"].index)))
only_T_1 = set(list(df_guias_nt[df_guias_nt["-1"] == "T"].index)).difference(set(list(df_guias_nt[df_guias_nt["+3"] == "A"].index)))

only_A_3 = set(list(df_guias_nt[df_guias_nt["+3"] == "A"].index)).difference(set(list(df_guias_nt[df_guias_nt["-1"] == "T"].index)))


venn2((len(only_T_1), len(only_A_3), len(inter)), set_labels = ('T in -1', 'A in +3'))
plt.title("Number of guides in each group")

plt.show()

inter_inserc = [guide_insertion[(x)] for x in list(inter)]
only_T_1_inserc = [guide_insertion[(x)] for x in list(only_T_1)]
only_A_3_inserc = [guide_insertion[(x)] for x in list(only_A_3)]

all_ins_data = inter_inserc + only_T_1_inserc + only_A_3_inserc

group = ["-1T and +3A"]*len(inter_inserc) +  ["Only -1T"]*len(only_T_1_inserc) +  ["Only +3A"]*len(only_A_3_inserc)
df_ins_by_group = pd.DataFrame({"Group":group, "Insertion":all_ins_data})

#sns.barplot(data = df_ins_by_group, x = "Group", y = "Insertion")
sns.boxplot(data = df_ins_by_group, x = "Group", y = "Insertion")
plt.title("Insertion percentage based on nucleotide ")
#sns.catplot(data = df_ins_by_group, x = "Group", y = "Insertion")

#%% All the possible combinations of 4 nucleotides (-2, -1, +1 y +3)
import itertools

seq = 'ACTG'
comb = []

#print(len(list(itertools.product(seq, seq, seq, seq))))

for output in itertools.product(seq, repeat=4):
    #print(''.join(output))
    comb.append(';'.join(output))
    

#%%

df_insertion_data_seq = pd.DataFrame(index = comb, columns = ["Mean insertion", "Number of guides", "Insertion", "Guides_ID"])

for k in comb: 
    nts = k.split(";")
    
    nt_l_2 = list(df_guias_nt[df_guias_nt["-2"] == nts[0]].index)
    nt_l_1 = list(df_guias_nt[df_guias_nt["-1"] == nts[1]].index)
    nt_r_1 = list(df_guias_nt[df_guias_nt["+1"] == nts[2]].index)
    nt_r_3 = list(df_guias_nt[df_guias_nt["+3"] == nts[3]].index)
    
    inter = set(nt_l_2) & set(nt_l_1) & set(nt_r_1) & set(nt_r_3)
    
    inter_inserc = [guide_insertion[(x)] for x in list(inter)]
    
    #print(k + " -> " + str(np.mean(inter_inserc)) + ". " + str(len(inter)) + " guías")
    df_insertion_data_seq.loc[k, :] = [np.median(inter_inserc), len(inter), inter_inserc, inter]

#%%

df_mayor_5 = df_insertion_data_seq[df_insertion_data_seq["Number of guides"] > 10]

df_mayor_5_sorted = (df_mayor_5.sort_values(by = "Mean insertion", ascending=False))
df_mayor_5_sorted.to_excel("..\\Insertion_data.xlsx")  
#%%

insertion = [item for sublist in list(df_mayor_5_sorted["Insertion"]) for item in sublist]
nucleotides = [[x]*len(df_mayor_5_sorted.loc[x, "Insertion"]) for x in list(df_mayor_5_sorted.index)]
nucleotides = [item for sublist in nucleotides for item in sublist]

df = pd.DataFrame(columns = ["Nucleotides", "Insertion"], index = nucleotides)
df["Nucleotides"] = nucleotides
df["Insertion"] = insertion

fig = sns.boxplot(data = df.loc[list(df_mayor_5_sorted.index), :], x = "Nucleotides", y = "Insertion")
fig.tick_params(labelsize=7)
fig.set_xticklabels(fig.get_xticklabels(),rotation=90)
plt.title("Distribution of insertion percentage \n in different nucleotide combinations")
#%%

fig = plt.figure(figsize=(10,6))
colors = ["y", "b", "g", "orange"]
c = 0

for nt in seq_list.keys():
    values = seq_list[nt]
    sns.distplot(values, color=colors[c],hist=False)
    c = c+1
 

fig.legend(labels=seq_list.keys(), loc = 'center right')
plt.title("Insertion percentage based on position -2", fontsize = 16)
plt.xlabel("% insertion")
plt.show()

#%%
from tabulate import tabulate
import random as rd

len(seq_list["G"])
len(seq_list["T"])
len(seq_list["A"])
len(seq_list["C"])



A_rn_data = rd.sample(seq_list["A"], k = 100)
T_rn_data = rd.sample(seq_list["T"], k = 100)
C_rn_data = rd.sample(seq_list["C"], k = 100)
G_rn_data = rd.sample(seq_list["G"], k = 100)


ins_mean = [['A', str(round(np.mean(seq_list["A"]),2)) + " ± " + str(round(np.std(seq_list["A"]),2))],
         ['C', str(round(np.mean(seq_list["C"]),2)) + " ± " + str(round(np.std(seq_list["C"]),2))],
         ['G', str(round(np.mean(seq_list["G"]),2)) + " ± " + str(round(np.std(seq_list["G"]),2))],
         ['T', str(round(np.mean(seq_list["T"]),2)) + " ± " + str(round(np.std(seq_list["T"]),2))]]
        
print(tabulate(ins_mean, headers=['Base en posición +3', '% Inserción']))


ins_mean_random = [['A', str(round(np.mean(A_rn_data),2)) + " ± " + str(round(np.std(A_rn_data),2))],
                   ['C', str(round(np.mean(C_rn_data),2)) + " ± " + str(round(np.std(C_rn_data),2))],
                   ['G', str(round(np.mean(G_rn_data),2)) + " ± " + str(round(np.std(G_rn_data),2))],
                   ['T', str(round(np.mean(T_rn_data),2)) + " ± " + str(round(np.std(T_rn_data),2))]]
                   
         

print(tabulate(ins_mean_random, headers=['Base en posición +3 \n (n = 100 random)', '% Inserción']))

"""
print("Porcentaje medio de inserción cuando hay una A en la posición -1: " + str(np.mean(seq_list["A"])))
print("Porcentaje medio de inserción cuando hay una C en la posición -1: " + str(np.mean(seq_list["C"])))
print("Porcentaje medio de inserción cuando hay una G en la posición -1: " + str(np.mean(seq_list["G"])))
print("Porcentaje medio de inserción cuando hay una T en la posición -1: " + str(np.mean(seq_list["T"])))
"""
#%%
#data = {'Mean_Ins': [np.mean(seq_list["A"]),np.mean(seq_list["T"]),np.mean(seq_list["C"]),np.mean(seq_list["G"])], 'Base': ["A", "T", "C", "G"]}

ins_data = [item for sublist in seq_list.values() for item in sublist]
base = [[x]*len(seq_list[x]) for x in seq_list.keys()]
base = [item for sublist in base for item in sublist]

df_mean_ins = pd.DataFrame({"Base":base, "Insertion":ins_data})

#sns.barplot(data = df_mean_ins, x = "Base", y = "Insertion")
sns.boxplot(data = df_mean_ins, x = "Base", y = "Insertion")
plt.title("Insertion percentage based on position +3")
#sns.catplot(data = df_mean_ins, x = "Base", y = "Insertion")

plt.show()
# Con los datos random

ins_data = A_rn_data + C_rn_data + T_rn_data + G_rn_data
base = ["A"]*100 +  ["C"]*100 +  ["T"]*100 +  ["G"]*100


df_mean_ins = pd.DataFrame({"Base":base, "Insertion":ins_data})

#sns.barplot(data = df_mean_ins, x = "Base", y = "Insertion")
sns.boxplot(data = df_mean_ins, x = "Base", y = "Insertion")
plt.title("Insertion percentage based on position +3")
plt.show()

#%% Entropy obtained with each target sequence

with open('..\\Logos\\file_sequence_entropy.txt', 'w') as f:
    for carpeta in ["results", "results_no_replic"]:
        
        os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta)
        guides_files = os.listdir()
        guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))
    
        for g_id in guides_id:
            
            #print(g_id)
            #print(glob.glob(str(g_id) + "_*"))
            files = glob.glob(str(g_id) + "_*")
            
            
            if len(files) >= 2:
                
                df_list_results = []
                df_list_insertion = []
    
                for file in files:  # si hay más de un archivo de la misma guía: 
                
                    df_list_results.append(pd.read_csv(file))
                    
                df_results = utils_functions.merge_df_results_alleles(df_list_results)
                df_results = df_results.fillna(0)
                
            else: 
                
                df_results = pd.read_csv(files[0])
                
            entropy_value = enp.entropy(df_results)
            sequence = (list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0]).upper()
            #sequence = (list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0])[24:34].upper()
            f.write(sequence + "\t" + str(entropy_value) + "\n")
            

#%%

df_entropy = pd.read_csv("..\\Logos\\file_sequence_entropy.txt", sep = "\t", names = ["Sequence", "Entropy"])

df_entropy["nt_3"] = [x[8] for x in df_entropy["Sequence"]]
#df_entropy["nt_3"] = [x[5] + x[8] for x in df_entropy["Sequence"]]

sns.boxplot(data = df_entropy.sort_values(by="Entropy", ascending = False), x = "nt_3",
            y = "Entropy")

plt.xlabel("Nucleotide at position +3")

plt.show()

len(df_entropy[df_entropy["nt_3"] == "T"])
len(df_entropy[df_entropy["nt_3"] == "G"])
len(df_entropy[df_entropy["nt_3"] == "A"])
len(df_entropy[df_entropy["nt_3"] == "C"])


import random as rd

A_rn_data = rd.sample(list(df_entropy[df_entropy["nt_3"] == "A"]["Entropy"]), k = 100)
T_rn_data = rd.sample(list(df_entropy[df_entropy["nt_3"] == "T"]["Entropy"]), k = 100)
C_rn_data = rd.sample(list(df_entropy[df_entropy["nt_3"] == "C"]["Entropy"]), k = 100)
G_rn_data = rd.sample(list(df_entropy[df_entropy["nt_3"] == "G"]["Entropy"]), k = 100)

df_entropy_rn100 = pd.DataFrame(columns =  ["Entropy","nt_3"])
df_entropy_rn100["Entropy"] = A_rn_data+T_rn_data+C_rn_data+G_rn_data
df_entropy_rn100["nt_3"] = ["A"]*100+["T"] * 100 + ["C"]*100 + ["G"] * 100

sns.boxplot(data = df_entropy_rn100.sort_values(by="Entropy", ascending = False), x = "nt_3", y = "Entropy")
plt.xlabel("Nucleotide at position +3")
plt.show()


from scipy import stats

a = df_entropy[df_entropy["nt_3"] == "A"]["Entropy"]
b = df_entropy[df_entropy["nt_3"] == "T"]["Entropy"]
t2, p2 = stats.ttest_ind(a,b)
print("t = " + str(t2))
print("p = " + str(p2))

a = list(df_entropy[df_entropy["nt_3"] == "T"]["Entropy"])
b = list(df_entropy[df_entropy["nt_3"] == "C"]["Entropy"])

a = T_rn_data
b = C_rn_data
