 # -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 16:38:48 2021

@author: yolib


Script of the general analysis of the alleles. 
We classify the 25 top alleles of each sample in different groups: 
    
    * D1 = 1bp deletions
    * D2 = 2bp deletions
    * D3 = 3bp deletions
    * D>3 = Deletions > 3bp
    * I1 = 1bp insertions
    * I2 = 2bp insertions
    * I>2 = Insertions > 2bp
    * Separated D = Two or more deletions in the same read
    * Separated I = Two or more insertions in the same read 
    * I+D = Insertion and deletions in the same read  
    
    
We examine too the size and the start position of the different indels. 
This analysis is performed using the guides with 3 or more replicates (526 guides)   
  
"""

#%%

import os
import pandas as pd
import numpy as np
import seaborn as sns 
import glob 
os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\python_scripts')
import utils_functions

#%% I load the data

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results')


guides_files = os.listdir()

lista_dict_min_freq = [] # diccionario de cada una de las muestras
lista_samples_min_freq = []
lista_guide_index_freq = []

for guide in guides_files: 
    data = pd.read_csv(guide)
    n_samples = data.shape[1] - 1
    #print(n_samples)
    #print(list(data.columns))
    #str(guide.split('.')[0] + )
    for sample in list(data.columns)[0:n_samples]:
        #print(guide.split('.')[0] + '_' + sample)
        sample_name_min_freq = str(guide.split('.')[0] + '_' + sample)
        sample_dict_min_freq = {}
        for k,line in data.iterrows():
            sample_dict_min_freq[line['Alleles']] = line[sample]
          
        lista_dict_min_freq.append(sample_dict_min_freq)    
        lista_samples_min_freq.append(sample_name_min_freq)
        lista_guide_index_freq.append(sample_name_min_freq.split('_')[0])



#%% Classification of the top 25 alleles. 


tipos_alelos = {}

I1 = []; I2 = []; I_higher_2 = []; I_sep = [] # I_sep son inserciones en lugares distintos
D1 = []; D2 = []; D3 = []; D_higher_3 = []; D_sep = []
I_D = []; SNV = []

#freq_min = 0.5 # voy a considerar solo aquellos que se encuentren en una frecuencia mayor del 0.5. 


for sample in lista_dict_min_freq:  
    
    #sample = lista_dict_min_freq[0]
    #alleles = set([x for x in sample if float(sample[x]) > freq_min])
    sample = dict(sorted(sample.items(), key=lambda item: item[1]))
    
    alleles = list(sample.keys())[-25:]
    
    if "no variant" in alleles: 
        alleles.remove('no variant')
        
    for allele in alleles: 
        
        if "SNV" in allele: 
            SNV.append(allele)
            
        else: 
            
            if "," in allele: 
                
                indels = [x[-1] for x in allele.split(",")]
                
                if "I" in indels and "D" in indels: 
                    I_D.append(allele)
                    
                elif "I" in indels and "D" not in indels: 
                    I_sep.append(allele)
                
                elif "D" in indels and "I" not in indels: 
                    D_sep.append(allele)
                    
            else: 
                
                if allele[-1] == "I": 
                    
                    if float(allele.split(":")[1][:-1]) == 1: 
                        I1.append(allele)
                    elif float(allele.split(":")[1][:-1]) == 2: 
                        I2.append(allele)
                    else: 
                        I_higher_2.append(allele)
                
                elif allele[-1] == "D": 
                    
                    if float(allele.split(":")[1][:-1]) == 1: 
                        D1.append(allele)
                    elif float(allele.split(":")[1][:-1]) == 2: 
                        D2.append(allele)
                    elif float(allele.split(":")[1][:-1]) == 3: 
                        D3.append(allele)
                    else: 
                        D_higher_3.append(allele)



tipos_alelos["D1"] = D1
tipos_alelos["D2"] = D2
tipos_alelos["D3"] = D3
tipos_alelos["D_higher_3"] = D_higher_3
tipos_alelos["I1"] = I1
tipos_alelos["I2"] = I2
tipos_alelos["I_higher_2"] = I_higher_2
tipos_alelos["D_sep"] = D_sep
tipos_alelos["I_sep"] = I_sep
tipos_alelos["I_D"] = I_D
tipos_alelos["SNV"] = SNV

#%%

total = 0

for tipo_al in tipos_alelos.keys():
    total = total + len(tipos_alelos[tipo_al])


#%%

dic_porc_tipo_alelos = {}

for tipo_al in tipos_alelos.keys():
    dic_porc_tipo_alelos[tipo_al] = 100*len(tipos_alelos[tipo_al])/total
    
#%% Pie chart of the top 25 alleles. 

import matplotlib.pyplot as plt
import random
import matplotlib.colors as mcolors
import seaborn as sns 

colors = sns.color_palette("hls", 14)

fig, ax = plt.subplots(figsize=(5, 6))

# Data to plot
labels = {}
porc = []

l = ["D1", "D2", "D3", "D>3", "I1", "I2", "I>2", "Separated D", "Separated I", "I+D", "SNV"]

c = 0
for x, y in dic_porc_tipo_alelos.items():
    labels[x] = l[c]
    porc.append(y)
    c = c+1


# Plot
plt.pie(porc, colors = colors,  radius = 1.0,
        #labels=labels, autopct='%1.1f%%',
        startangle=90)

ax.legend([str(labels[x] + ": " + str(round(dic_porc_tipo_alelos[x], 2)) + "%") for x in dic_porc_tipo_alelos.keys()],
          title="Alleles",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))

#plt.title("Proporci??n de los distintos alelos con frecuencia \n mayor a " + str(freq_min) + "%")
#plt.title("Proportion of the alleles with a frequency \n higher than " + str(freq_min) + "%")
plt.title("Proportion of the alleles found among the top 25")
#plt.axis('equal')
plt.show()


#%% Insertions' size


fig, ax = plt.subplots(2)

size_1 = [1] * len(tipos_alelos["I1"])
size_2 = [2] * len(tipos_alelos["I2"])
size_mas_2 = [int(x.split(':')[1][:-1]) for x in tipos_alelos["I_higher_2"]]

size_I = size_1 + size_2 + size_mas_2

sns.histplot(size_I,discrete=True, ax = ax[0], stat = "probability",color="cornflowerblue")
#plt.title("Frecuencia del tama??o de inserciones")
plt.xlabel("Indel size", size = 2)
#ax[0].set_xticks(np.arange(0, max(size_I), 5))
#ax[0].set(xlabel='tama??o inserciones')
          #, ylabel='y-label')
fig.suptitle("Size distribution of indels", fontsize = 15)


size_I = np.array(size_I)
print("Probabilidad de inserci??n de 1 base: " + str(len(size_I[size_I == 1])/len(size_I)))

# An??lisis del tama??o de las deleciones

size_D1 = [1] * len(tipos_alelos["D1"])
size_D2 = [2] * len(tipos_alelos["D2"])
size_D3 = [3] * len(tipos_alelos["D3"])
size_mas_D3 = [int(x.split(':')[1][:-1]) for x in tipos_alelos["D_higher_3"] ] #if (int(x.split(':')[1][:-1]) + int(x.split(':')[0])) < 60]

size_D = size_D1 + size_D2 + size_D3 + size_mas_D3

sns.histplot(size_D,discrete=True, ax = ax[1], stat = "probability",color="cornflowerblue")
ax[1].set_xticks(np.arange(0, max(size_D), 5))
#ax[1].set(xlabel='tama??o deleciones')
plt.xlabel("Indel size", size = 12)


ax[0].set_title('Insertions', fontdict = {'fontsize' : 13})
ax[1].set_title('Deletions', fontdict = {'fontsize' : 13})

plt.tight_layout()

#%% Starting positions of the insertions

fig, ax = plt.subplots(2)

pos_1 = [int(x.split(':')[0]) for x in tipos_alelos["I1"]]
pos_2 = [int(x.split(':')[0]) for x in tipos_alelos["I2"]]
pos_more_2 = [int(x.split(':')[0]) for x in tipos_alelos["I_higher_2"]]

pos_I = pos_1 + pos_2 + pos_more_2


sns.histplot(pos_I,discrete=True, ax = ax[0], stat = "probability",color="cornflowerblue")
#plt.title("Frecuencia del tama??o de inserciones")
plt.xlabel("Indel position", size = 2)
#ax[0].set_xticks(np.arange(0, max(size_I), 5))
#ax[0].set(xlabel='tama??o inserciones')
          #, ylabel='y-label')
fig.suptitle("Distribution of the starting position of the indels", fontsize = 15)

pos_I = np.array(pos_I)
a = len(pos_I[pos_I == 1])
b = len(pos_I[pos_I == 2])
c = len(pos_I[pos_I == 3])
d = len(pos_I[pos_I == -1])
e = len(pos_I[pos_I == -2])

(a+b+c+d+e)/len(pos_I)

print(a+b+c+d+e)

# An??lisis de las posiciones de las deleciones

posD1 = [int(x.split(':')[0]) for x in tipos_alelos["D1"]]
posD2 = [int(x.split(':')[0]) for x in tipos_alelos["D2"]]
posD3 = [int(x.split(':')[0]) for x in tipos_alelos["D3"]]
pos_more_D3 = [int(x.split(':')[0]) for x in tipos_alelos["D_higher_3"]]

pos_D = posD1 + posD2 + posD3 + pos_more_D3

sns.histplot(pos_D,discrete=True, ax = ax[1], stat = "probability",color="cornflowerblue")
#ax[1].set_xticks(np.arange(0, max(size_D), 5))
#ax[1].set(xlabel='tama??o deleciones')
plt.xlabel("Indel position", size = 12)


ax[0].set_title('Insertions', fontdict = {'fontsize' : 13})
ax[1].set_title('Deletions', fontdict = {'fontsize' : 13})

plt.tight_layout()


#%% Insertions of 1 bp


pos_I1 = [int(x.split(':')[0]) for x in tipos_alelos["I1"]]

#plt.hist(pos_I1, bins = np.arange(min(pos_I1), max(pos_I1), 1), alpha=1, edgecolor = 'black',  linewidth=1)
fig_dims = (8, 6)
fig, ax = plt.subplots(figsize=fig_dims)

sns.histplot(pos_I1,discrete=True, stat = "probability", ax = ax, color="cornflowerblue")
plt.title("Position of 1 base inserts", fontsize = 15)
plt.xlabel("Position with respect to the cleavage site ", size = 12)
"""
plt.title("Posici??n de las inserciones de 1 base")
plt.xlabel("Posici??n con respecto al punto de corte", size = 12)
"""
#%% Insertions of 2 bp

pos_I2 = [int(x.split(':')[0]) for x in tipos_alelos["I2"]]

fig_dims = (8, 6)
fig, ax = plt.subplots(figsize=fig_dims)

sns.histplot(pos_I2,discrete=True, stat = "probability", ax = ax, color="cornflowerblue")
#plt.title("Position of 2 bp insertions", fontsize = 15)
plt.xlabel("Position with respect to the cleavage site ", size = 12)

#%% Deletions of 1 bp

pos_D1 = [int(x.split(':')[0]) for x in tipos_alelos["D1"]]

fig_dims = (8, 6)
fig, ax = plt.subplots(figsize=fig_dims)

sns.histplot(pos_D1,discrete=True, stat = "probability", ax = ax, color="cornflowerblue")
plt.title("Position of 1 bp deletions", fontsize = 15)
plt.xlabel("Position with respect to the cleavage site ", size = 12)


#%% Deletions of 2 bp 

pos_D2 = [int(x.split(':')[0]) for x in tipos_alelos["D2"]]

fig_dims = (8, 6)
fig, ax = plt.subplots(figsize=fig_dims)

sns.histplot(pos_D2,discrete=True, stat = "probability", ax = ax, color="cornflowerblue")
#plt.title("Position of 2 bp deletions", fontsize = 15)
plt.xlabel("Position with respect to the cleavage site ", size = 12)

#%% Deletions of 3 bp

fig_dims = (8, 6)
fig, ax = plt.subplots(figsize=fig_dims)

pos_D3 = [int(x.split(':')[0]) for x in tipos_alelos["D3"]]

sns.histplot(pos_D3,discrete=True, stat = "probability", ax = ax, color="cornflowerblue")
#plt.title("Position of 3 bp deletions")
plt.xlabel("Position with respect to the cleavage site ", size = 12)
ax.set_xticks(np.arange(-60, max(pos_D3), 5))
plt.xticks(fontsize=10)

#%% Deletions of > 3 bp

fig_dims = (8, 6)
fig, ax = plt.subplots(figsize=fig_dims)

pos_D_mayor_3 = [int(x.split(':')[0]) for x in tipos_alelos["D_higher_3"]]

sns.histplot(pos_D_mayor_3,discrete=True, stat = "probability", ax = ax, color="cornflowerblue")
plt.title("Position of deletions larger than 3 bp", size = 12)
plt.xlabel("Position with respect to the cleavage site ", size = 12)
ax.set_xticks(np.arange(-65, max(pos_D_mayor_3), 5))
plt.xticks(fontsize=8)
#%% Sizes of the deletions > 3bp

fig_dims = (8, 6)
fig, ax = plt.subplots(figsize=fig_dims)

size_D_mayor_3 = [int(x.split(':')[1][:-1]) for x in tipos_alelos["D_higher_3"]]

sns.histplot(size_D_mayor_3,discrete=True, ax = ax , color="cornflowerblue",stat = "probability",)
plt.title("Size of deletions larger than 3 bp",size = 12,)
plt.xlabel("Deletion size", size = 12)

ax.set_xticks(np.arange(0, max(size_D_mayor_3), 5))
plt.xticks(fontsize=10) 


#%% Most common allele

import operator

dict_alelo_mayorit = {}

I1 = []; I2 = []; I_mayor_2 = []; I_sep = [] # I_sep son inserciones en lugares distintos
D1 = []; D2 = []; D3 = []; D_mayor_3 = []; D_sep = []
I_D = []; SNV = []

for sample in lista_dict_min_freq: 
    
    sorted_alleles = sorted(sample.items(), key=operator.itemgetter(1), reverse=True)    
    
    if "no variant" in sorted_alleles[0]: 
        alelo_mayor = sorted_alleles[1][0]
    else: 
        alelo_mayor = sorted_alleles[0][0]
        
    if "SNV" in alelo_mayor: 
        SNV.append(alelo_mayor)
        
    elif "," in alelo_mayor: 
        
        indels = [x[-1] for x in alelo_mayor.split(",")]
        
        if "I" in indels and "D" in indels: 
            I_D.append(alelo_mayor)
            
        elif "I" in indels and "D" not in indels: 
            I_sep.append(alelo_mayor)
        
        elif "D" in indels and "I" not in indels: 
            D_sep.append(alelo_mayor)
            
    else: 
        
        if alelo_mayor[-1] == "I": 
            
            if float(alelo_mayor.split(":")[1][:-1]) == 1: 
                I1.append(alelo_mayor)
            elif float(alelo_mayor.split(":")[1][:-1]) == 2: 
                I2.append(alelo_mayor)
            else: 
                I_mayor_2.append(alelo_mayor)
        
        elif alelo_mayor[-1] == "D": 
            
            if float(alelo_mayor.split(":")[1][:-1]) == 1: 
                D1.append(alelo_mayor)
            elif float(alelo_mayor.split(":")[1][:-1]) == 2: 
                D2.append(alelo_mayor)
            elif float(alelo_mayor.split(":")[1][:-1]) == 3: 
                D3.append(alelo_mayor)
            else: 
                D_mayor_3.append(alelo_mayor)



dict_alelo_mayorit["D1"] = D1
dict_alelo_mayorit["D2"] = D2
dict_alelo_mayorit["D3"] = D3
dict_alelo_mayorit["D_mayor_3"] = D_mayor_3
dict_alelo_mayorit["I1"] = I1
dict_alelo_mayorit["I2"] = I2
dict_alelo_mayorit["I_mayor_2"] = I_mayor_2
dict_alelo_mayorit["D_sep"] = D_sep
dict_alelo_mayorit["I_sep"] = I_sep
dict_alelo_mayorit["I_D"] = I_D
dict_alelo_mayorit["SNV"] = SNV

#%%
total_mayor = 0

for mayor_al in dict_alelo_mayorit.keys():
    total_mayor = total_mayor + len(dict_alelo_mayorit[mayor_al])

#%%

dic_porc_mayor_alelo = {}

for mayor_al in dict_alelo_mayorit.keys():
    dic_porc_mayor_alelo[mayor_al] = 100*len(dict_alelo_mayorit[mayor_al])/total_mayor


#%%

import matplotlib.pyplot as plt
import random
import matplotlib.colors as mcolors
import seaborn as sns 

colors = sns.color_palette("hls", 14)

fig, ax = plt.subplots(figsize=(5, 6))


# Data to plot
labels = {}
porc = []

l = ["D1", "D2", "D3", "D>3", "I1", "I2", "I>2", "Separated D", "Separated I", "I+D", "SNV"]

c = 0
for x, y in dic_porc_mayor_alelo.items():
    labels[x] = l[c]
    porc.append(y)
    c = c+1


# Plot
plt.pie(porc, colors = colors, radius = 1.0,
        #labels=labels, #autopct='%1.1f%%',
        startangle=90)

ax.legend([str(labels[x] + ": " + str(round(dic_porc_mayor_alelo[x], 2)) + "%") for x in dic_porc_mayor_alelo.keys()],
          title="Alleles",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))

#plt.title("Proporci??n de los distintos alelos con frecuencia \n mayor a " + str(freq_min) + "%")
#plt.title("Frequency of the most common allele  ")

plt.axis('equal')
plt.show()

# This is all with respect to the majority allele
#%% Insertions of 1bp

pos_I1 = [int(x.split(':')[0]) for x in dict_alelo_mayorit["I1"]]

#plt.hist(pos_I1, bins = np.arange(min(pos_I1), max(pos_I1), 1), alpha=1, edgecolor = 'black',  linewidth=1)

sns.displot(pos_I1,discrete=True)
plt.title("Posici??n de las inserciones de 1 base (alelo mayoritario)")
plt.xlabel("Posici??n con respecto al punto de corte", size = 12)


#%% Insertions of 2bp
"""
pos_I2 = [int(x.split(':')[0]) for x in dict_alelo_mayorit["I2"]]

sns.displot(pos_I2,discrete=True)
plt.title("Frecuencia de inserci??n I2")
plt.xlabel("Posici??n con respecto al punto de corte", size = 8)
"""
#%% 1bp deletions

pos_D1 = [int(x.split(':')[0]) for x in dict_alelo_mayorit["D1"]]

sns.displot(pos_D1,discrete=True)
plt.title("Posici??n de las deleciones de 1 base (alelo mayoritario)")
plt.xlabel("Posici??n con respecto al punto de corte", size = 12)
#%% 2bp deletions

pos_D2 = [int(x.split(':')[0]) for x in dict_alelo_mayorit["D2"]]

sns.displot(pos_D2,discrete=True)
plt.title("Posici??n de las deleciones de 2 bases (alelo mayoritario)")
plt.xlabel("Posici??n con respecto al punto de corte", size = 12)

#%%

#%% Most common allele replicable between 3 or more replicas


from collections import Counter

mayor_allele_replicas = {} 


I1 = []; I2 = []; I_mayor_2 = []; I_sep = [] # I_sep son inserciones en lugares distintos
D1 = []; D2 = []; D3 = []; D_mayor_3 = []; D_sep = []
I_D = []; SNV = []
no_common = []; more_than_one_common = []

guides = list(set(lista_guide_index_freq))



for guide in guides: 
    
    lista_al_mayor = []
    
    indices_guides = [i for i,d in enumerate(lista_guide_index_freq) if d==guide]
    
    for sample_i in indices_guides: 
        
        sorted_alleles = sorted(lista_dict_min_freq[sample_i].items(), key=operator.itemgetter(1), reverse=True) 
        
        if "no variant" in sorted_alleles[0]: 
            alelo_mayor = sorted_alleles[1][0]
            #corte = sorted_alleles[1][1]
        else: 
            alelo_mayor = sorted_alleles[0][0]
            #corte = sorted_alleles[0][1]
    
        lista_al_mayor.append(alelo_mayor)

    
    rep = Counter(lista_al_mayor)
    
    if len(list(rep)) == 1: 
        
        if "SNV" in list(rep)[0]: 
            SNV.append(list(rep)[0])
        
        if ',' in list(rep)[0]: 
            
            indels = [x[-1] for x in list(rep)[0].split(",")]
        
            if "I" in indels and "D" in indels: 
                I_D.append(list(rep)[0])
                
            elif "I" in indels and "D" not in indels: 
                I_sep.append(list(rep)[0])
            
            elif "D" in indels and "I" not in indels: 
                D_sep.append(list(rep)[0])
            
        else: 
            
            if list(rep)[0][-1] == "I": 
            
                if float(list(rep)[0][-2]) == 1: 
                    I1.append(list(rep)[0])
                elif float(list(rep)[0][-2]) == 2: 
                    I2.append(list(rep)[0])
                else: 
                    I_mayor_2.append(list(rep)[0])
        
            elif list(rep)[0][-1] == "D": 
                
                if float(list(rep)[0][-2]) == 1: 
                    D1.append(list(rep)[0])
                elif float(list(rep)[0][-2]) == 2: 
                    D2.append(list(rep)[0])
                elif float(list(rep)[0][-2]) == 3: 
                    D3.append(list(rep)[0])
                else: 
                    D_mayor_3.append(list(rep)[0])
            
            
    else: # si no todas las r??plicas tienen el mismo alelo mayoritario: 
        rep = rep.most_common()

        if rep[0][1] > 3 and rep[0][1] != rep[1][1]: 
            
            if ',' in list(rep)[0][0]: 
            
                indels = [x[-1] for x in list(rep)[0][0].split(",")]
            
                if "I" in indels and "D" in indels: 
                    I_D.append(list(rep)[0][0])
                    
                elif "I" in indels and "D" not in indels: 
                    I_sep.append(list(rep)[0][0])
                
                elif "D" in indels and "I" not in indels: 
                    D_sep.append(list(rep)[0][0])
            
            else: 
                
                if list(rep)[0][0][-1] == "I": 
                
                    if float(list(rep)[0][0][-2]) == 1: 
                        I1.append(list(rep)[0][0])
                    elif float(list(rep)[0][0][-2]) == 2: 
                        I2.append(list(rep)[0][0])
                    else: 
                        I_mayor_2.append(list(rep)[0][0])
        
                elif list(rep)[0][0][-1] == "D": 
                    
                    if float(list(rep)[0][0][-2]) == 1: 
                        D1.append(list(rep)[0][0])
                    elif float(list(rep)[0][0][-2]) == 2: 
                        D2.append(list(rep)[0][0])
                    elif float(list(rep)[0][0][-2]) == 3: 
                        D3.append(list(rep)[0][0])
                    else: 
                        D_mayor_3.append(list(rep)[0][0])
        
                elif "SNV" in list(rep)[0][0]: 
                    SNV.append(list(rep)[0][0])
        
        elif rep[0][1] >= 3 and rep[0][1] == rep[1][1]: 
    
            more_than_one_common.append(rep[0][0])

        elif rep[0][1] < 3 : 
            
            no_common.append(rep)


mayor_allele_replicas["D1"] = D1
mayor_allele_replicas["D2"] = D2
mayor_allele_replicas["D3"] = D3
mayor_allele_replicas["D_mayor_3"] = D_mayor_3
mayor_allele_replicas["I1"] = I1
mayor_allele_replicas["I2"] = I2
mayor_allele_replicas["I_mayor_2"] = I_mayor_2
mayor_allele_replicas["D_sep"] = D_sep
mayor_allele_replicas["I_sep"] = I_sep
mayor_allele_replicas["I_D"] = I_D
mayor_allele_replicas["SNV"] = SNV
mayor_allele_replicas["no_comun"] = no_common # When at least 3 replicates don't have the same most common allele
mayor_allele_replicas["mas_1_comun"] = more_than_one_common # If there is more than 1 most common allele (for example, 3 have one and other 3 another)

#%%
total_mayor = 0

for mayor_al in mayor_allele_replicas.keys():
    total_mayor = total_mayor + len(mayor_allele_replicas[mayor_al])

#%%

dic_porc_mayor_alelo_comun = {}

for mayor_al in mayor_allele_replicas.keys():
    dic_porc_mayor_alelo_comun[mayor_al] = 100*len(mayor_allele_replicas[mayor_al])/total_mayor


#%% 
import matplotlib.pyplot as plt
import random
import matplotlib.colors as mcolors
import seaborn as sns 

colors = sns.color_palette("hls", 14)

fig, ax = plt.subplots(figsize=(5, 6))


# Data to plot
labels = {}
porc = []

l = ["D1", "D2", "D3", "D>3", "I1", "I2", "I>2", "Separated D", "Separated I", "I+D", "SNV", "No common", "More than 1 common"]

c = 0
for x, y in dic_porc_mayor_alelo_comun.items():
    labels[x] = l[c]
    porc.append(y)
    c = c+1


# Plot
plt.pie(porc, colors = colors, radius= 1.0,
        #labels=labels, #autopct='%1.1f%%',
        startangle=90)

ax.legend([str(labels[x] + ": " + str(round(dic_porc_mayor_alelo_comun[x], 2)) + "%") for x in dic_porc_mayor_alelo_comun.keys()],
          title="Alleles",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))

#plt.title("Proporci??n de los distintos alelos con frecuencia \n mayor a " + str(freq_min) + "%")
#plt.title("Frequency of the common majority allele \n among 3 or more replicates")

plt.axis('equal')
plt.show()

