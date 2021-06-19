# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 17:43:10 2021

@author: yolib

Function to merge the dataframes of experiments resulting from the same guide. 
"""


def merge_df_results_alleles(df_list): 
    # function to join the dataframes according to the alleles column. 
    # If an allele is not in the other replicas, it writes NaN.
    
    import numpy as np
    
    if len(df_list) == 2: 
        
        merged = df_list[0].merge(df_list[1],how='outer', on = "Alleles")
        
    elif len(df_list) > 2: 
        
        merged = df_list[0].merge(df_list[1],how='outer', on = "Alleles")
            
        for df_i in np.arange(2,len(df_list)): 
            
            merged = merged.merge(df_list[df_i],how='outer', on = "Alleles")
        
    
    col_names = list(merged.columns)
    col_names.remove("Alleles")
    merged = merged[col_names + ["Alleles"]]
    
    return merged