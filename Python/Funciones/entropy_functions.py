# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 13:08:32 2021

@author: yolib

Function to compute the Shannon entropy. 

"""

def entropy(df_results): 
    
    # Function to compute the Shannon Entropy of a mutational profile. 
    # The input is the dataframe of the mutational profile. 
    
    import numpy as np 
    import math 
    
    if "no variant" in list(df_results["Alleles"]): 
        porc_corte = 100-np.mean(list(df_results[df_results["Alleles"] == "no variant"].iloc[0, :-1]))
        #df_results = df_results.drop(df_results[df_results['Alleles']=="no variant"])
        df_results = df_results[df_results.Alleles != 'no variant']
    else: 
        porc_corte = 100
        
    entropy_value = 0
    
    for row in list(df_results.index):
        #print(row)
        if "," not in df_results.loc[row, "Alleles"] and "SNV" not in df_results.loc[row, "Alleles"]: 
            pi = np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]]))/porc_corte
            entropy_value = entropy_value + (pi * math.log(pi,2))

    return entropy_value*(-1)        


def entropy_25(df_results): 
    # Function to compute the Shannon Entropy of a mutational profile. 
    # The input is the dataframe of the mutational profile. 
    # This function only takes the 25 top alleles (the 25 most frequent alleles). 
    
    import numpy as np 
    import math 
    
    if "no variant" in list(df_results["Alleles"]): 
        porc_corte = 100-np.mean(list(df_results[df_results["Alleles"] == "no variant"].iloc[0, :-1]))
        #df_results = df_results.drop(df_results[df_results['Alleles']=="no variant"])
        df_results = df_results[df_results.Alleles != 'no variant']
    else: 
        porc_corte = 100
        
    entropy_value = 0
    df_results = df_results.sort_values(by = list(df_results.columns)[0:-1], ascending = False)
                            
    for row in list(df_results.index)[0:25]:
        #print(row)
        pi = np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]]))/porc_corte
        entropy_value = entropy_value + (pi * math.log(pi,2))

    return entropy_value*(-1)      


