# -*- coding: utf-8 -*-
"""
Created on Mon May 17 13:38:50 2021

@author: yolib

Script to compute the MH score of all the guides from the dataset of Leenay et al. 

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
os.chdir('C:\\Users\\yolib\\OneDrive\\Documentos\\TFM\\Documento TFM\\Codigo\\Python\\Funciones')
import MH_functions 
import utils_functions

#%%


os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos')
guides_data_no_replic = pd.read_csv("guides_sequence_2_no_replic.csv")
guides_replic_data = pd.read_csv("guides_sequence_2.csv")

frames = [guides_replic_data, guides_data_no_replic]
guides_data = pd.concat(frames, ignore_index = True)

guides_data = guides_data.drop_duplicates()


#%% We compute the pattern scores, the MH scores and the out of frame scores. 

MH_score_dict = {}
out_of_frame_score_dict = {}
pattern_score_most_prob = {}

for guide_ID in list(guides_data["Guide_ID"]): 
    sequence = list(guides_data[guides_data["Guide_ID"] == guide_ID]["Sequence"])[0]
    MH_functions.MH_patterns(sequence.upper(), 29, file = 'C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\MH\\' + str(guide_ID) + ".txt")
    df_pattern_score = MH_functions.pattern_score(file = 'C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\MH\\' + str(guide_ID) + ".txt")
    pattern_score_most_prob[guide_ID] = df_pattern_score.iloc[0, -1]
    MH_score_dict[guide_ID] = df_pattern_score["Pattern_Score"].sum()
    out_of_frame_score_dict[guide_ID] = 100*df_pattern_score[df_pattern_score["Out_of_frame"] == True]["Pattern_Score"].sum()/MH_score_dict[guide_ID]
    
    df_pattern_score.to_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\MH\\' + str(guide_ID) + ".txt", sep = "\t")


#%% We take the 25 top alleles and we calculate the percentage of MH deletions and the percentage of out of frame deletions. 

porc_MH_deletion_dict = {}
porc_out_frame_dict = {}

for carpeta in ["results", "results_no_replic"]:
    
    os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta)
    guides_files = os.listdir()
    guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))

    for g_id in guides_id:
        
        #print(g_id)
        #print(glob.glob(str(g_id) + "_*"))
        files = glob.glob(str(g_id) + "_*")
        
       
        #list_MH_del = []
        
        if len(files) >= 2:
            
            df_list_results = []
            df_list_insertion = []

            for file in files:  # si hay más de un archivo de la misma guía: 
            
                df_list_results.append(pd.read_csv(file))
                
            df_results = utils_functions.merge_df_results_alleles(df_list_results)
            df_results = df_results.fillna(0)
            
        else: 
            
            df_results = pd.read_csv(files[0])
            
        if "no variant" in list(df_results["Alleles"]):
            porc_corte = (100-np.mean(list(df_results[df_results["Alleles"] == "no variant"].iloc[0, :-1])))
            #porc_corte_list.append(porc_corte)
            df_results = df_results[df_results["Alleles"] != "no variant"]
        
        else: 
            porc_corte = 100
            
        norm = 1/porc_corte
        
        if porc_corte >= 10: # porque si el % de corte es muy pequeño, me mete mucho ruido
            
            sequence = list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0]
            
            df_results = df_results.sort_values(by = list(df_results.columns)[:-1], ascending = False)
            #df_results = df_results.iloc[0:25, :]    
            
            alleles_deletion = [x for x in list(df_results["Alleles"]) if "," not in x and x[-1] == "D"]
            deletions_index = list(df_results.loc[df_results['Alleles'].isin(alleles_deletion)].index)
            
            porc_deletion = 0.0
            
            for row in deletions_index:
                #if np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]])) >= 0.5:
                porc_deletion = porc_deletion + np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]]))
           
            
            out_of_frame_deletions = [x for x in alleles_deletion if int(x.split(":")[1][:-1])%3 != 0]
            out_of_frame_deletions_index = list(df_results.loc[df_results['Alleles'].isin(out_of_frame_deletions)].index)
            
            porc_out_frame = 0.0
            
            for row in out_of_frame_deletions_index:
                #if np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]])) >= 0.5: 
                porc_out_frame = porc_out_frame + np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]]))
           
            
            list_MH_del = MH_functions.all_MH_cigar(sequence.upper(), alleles_deletion)
            MH_deletions_index = list(df_results.loc[df_results['Alleles'].isin(list_MH_del)].index)
           
            porc_MH_deletion = 0.0
            
            for row in MH_deletions_index:
                #if np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]])) >= 0.5: 
                #porc_MH_deletion = porc_MH_deletion + np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]]))
                
                porc_MH_deletion = porc_MH_deletion + np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]]))
            
            if porc_deletion > 0:
                porc_MH_deletion_dict[g_id] = porc_MH_deletion/porc_deletion
                porc_out_frame_dict[g_id] = porc_out_frame/porc_deletion


#%%

df_data = pd.DataFrame(index = porc_MH_deletion_dict.keys(), 
                                 columns = ["MH_score", "MH_deletion"])
            
for k in porc_MH_deletion_dict.keys(): 
    df_data.loc[k, "MH_score"] = float(MH_score_dict[k])  
    df_data.loc[k, "MH_deletion"] = float(porc_MH_deletion_dict[k]  )



fig, ax = plt.subplots(figsize=(15, 10))

sns.scatterplot(x = "MH_score",
                y = "MH_deletion", data = df_data,
                ax = ax)


#%% Correlation between the MH score and the % of deletions mediated by MH. 


from numpy.random import randn
from numpy.random import seed
from scipy.stats import pearsonr
import scipy.stats 
from matplotlib.font_manager import FontProperties

x = list(df_data["MH_score"])
y = list(df_data["MH_deletion"])

slope, intercept, r, p, stderr = scipy.stats.linregress(x, y)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fontP = FontProperties()
fontP.set_size('medium')
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(x, y, linewidth=0, marker='o', label='Data points', 
        markersize=5, c='cornflowerblue')
ax.plot(x, intercept + slope * np.array(x), 
        label=line, c = "black", linestyle=":")
ax.set_ylabel("% MH deletions", size = 12)
ax.set_xlabel("MH score", size = 12)
plt.title("Correlation between the MH score \n and the proportion of deletions produced by MH", size = 15)
ax.legend(facecolor='white', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
text = str('R = ') + str(round(r,3))
ax.text(1.2, 0.85, text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.show()


#%% 

df_data_out_frame = pd.DataFrame(index = porc_out_frame_dict.keys(), 
                                 columns = ["MH_out_frame_score", "Porc_out_frame"])
            
for k in porc_out_frame_dict.keys(): 
    df_data_out_frame.loc[k, "MH_out_frame_score"] = float(out_of_frame_score_dict[k])  
    df_data_out_frame.loc[k, "Porc_out_frame"] = float(porc_out_frame_dict[k]  )



fig, ax = plt.subplots(figsize=(15, 10))

sns.scatterplot(x = "MH_out_frame_score",
                y = "Porc_out_frame", data = df_data_out_frame,
                ax = ax)


#%% Correlation between the out of frame score and the % of deletions out of frame. 



from numpy.random import randn
from numpy.random import seed
from scipy.stats import pearsonr
import scipy.stats 
from matplotlib.font_manager import FontProperties

y = list(df_data_out_frame["Porc_out_frame"])
x = list(df_data_out_frame["MH_out_frame_score"])

slope, intercept, r, p, stderr = scipy.stats.linregress(x, y)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fontP = FontProperties()
fontP.set_size('medium')
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(x, y, linewidth=0, marker='o', label='Data points', 
        markersize=5, c='cornflowerblue')
ax.plot(x, intercept + slope * np.array(x), 
        label=line, c = "black", linestyle=":")
ax.set_xlabel("Out of frame score", size = 12)
ax.set_ylabel("Proportion of out of frame deletions", size = 12)
plt.title("Correlation between the MH out of frame score \n and the proportion of out of frame deletions", size = 15)
ax.legend(facecolor='white', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
text = str('R = ') + str(round(r,3))
ax.text(1.2, 0.85, text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.show()

