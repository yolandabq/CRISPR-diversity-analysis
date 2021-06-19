# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 18:10:09 2021

@author: yolib

Functions to analyze if the repairing process is reproducible.
Functions to calculate the symmetric Kullback-Leibler divergence and the % of in frame mutations. 

"""


def KL(p1, p2, missing = 0.001): 
    
    """
    Function to calculate the divergence of KL of one sample with respect to another. This is the
    asymmetric divergence.
    
    p1 is the pattern of mutations from sample 1 and p2 is from sample 2. p1 and p2 are dictionaries
    whose key are the different alleles and their values ​​are their corresponding frequencies.
    
    "missing" is the value of the pseudo-account that we want to add to avoid dividing by 0 when an allele is found in
    one sample but not another (since its frequency would be 0).
    
    """
    
    from math import log2
    
    p1_alleles = set([x for x in p1 if p1[x]>0])
    p2_alleles = set([x for x in p2 if p2[x]>0])
    
    if "no variant" in p1_alleles: 
        p1_alleles.remove('no variant')
        
    if "no variant" in p2_alleles: 
        p2_alleles.remove('no variant')
        
    
    common_alleles = p1_alleles.intersection(p2_alleles)
    p1_only_alleles = p1_alleles.difference(p2_alleles)
    p2_only_alleles = p2_alleles.difference(p1_alleles)
    
    p1_suma = sum([p1[allele] for allele in p1_alleles]) + missing*len(p2_only_alleles)
    p2_suma = sum([p2[allele] for allele in p2_alleles]) + missing*len(p1_only_alleles)
    
    if "no variant" in p1.keys(): 
        p1_porc_corte = 100 - p1['no variant']
        p1_norm = 1/p1_porc_corte # normalizamos dividiendo entre el % de corte. 
    else: 
        p1_norm = 1/p1_suma
        
    if "no variant" in p2.keys(): 
        p2_porc_corte = 100 - p2['no variant']
        p2_norm = 1/p2_porc_corte
    else: 
        p2_norm = 1/p2_suma
        
    KL_score = 0
    
    for mut in common_alleles: 
        KL_score += p1[mut]*p1_norm*log2(p1[mut]*p1_norm/(p2[mut]*p2_norm))
    for mut in p1_only_alleles: 
         KL_score += p1[mut]*p1_norm*log2(p1[mut]*p1_norm/(missing*p2_norm))
    for mut in p2_only_alleles: 
        KL_score += missing*p1_norm*log2(missing*p1_norm/(p2[mut]*p2_norm))
    
    return KL_score


def symmetricKL(p1, p2,  missing = 0.001): # function to calculate the KL of the two samples (p1 versus p2 and p2 versus p1,
# which gives different results, and we add them together)
    """
    
    This function calculates the symmetric KL divergence. To do this, it obtains the value of the asymmetric divergence of KL
    of p1 over p2 and of p2 over p1 and add both values.
    
    p1 is the pattern of mutations from sample 1 and p2 is from sample 2. p1 and p2 are dictionaries
    whose key are the different alleles and their values ​​are their corresponding frequencies.
    
    "missing" is the value of the pseudo-account that we want to add to avoid dividing by 0 when an allele is found in
    one sample but not another (since its frequency would be 0).
    
    """
    return 0.5*KL(p1, p2, missing) + 0.5*KL(p2, p1, missing)
    

def in_frame(sample):
    
    #print('---')
    """
    Function to determine the percentage of mutations "in frame", that is, that
    do not lose the reading frame (they have a size of 3 or multiple of 3)
    The input is a dictionary with the different mutations (keys) and their percentage
    (value of each key). For example: sample [1:1I] = 30
    """
    
    alleles = list(set([x for x in sample]))
    alleles_in_frame = []
    frame_shift_alleles = []
    
    if "no variant" in alleles: 
        porc_corte = 100 - sample["no variant"]
        alleles.remove('no variant')
        
        
    else: 
        porc_corte = 100
        
    norm = 1/porc_corte
    
    for a in alleles: 
        
        if "SNV" in a: 
            alleles_in_frame.append(sample[a] * norm)
            #pass
        
        elif "," in a: 
            size = 0
            
            for k in a.split(","): 
                
                size = size + int(k.split(":")[1][:-1])
                
            if size % 3 != 0: # If it isn't multiple of 3
                
                frame_shift_alleles.append(sample[a] * norm)
                
            else: 
                
                alleles_in_frame.append(sample[a] * norm)
    
        else: 
            
            if int(a.split(":")[1][:-1]) % 3 != 0: 
                
                frame_shift_alleles.append(sample[a] * norm)
                
            else: 
                
                alleles_in_frame.append(sample[a] * norm)
    
    in_frame_porc = sum(alleles_in_frame)
    out_frame_porc = sum(frame_shift_alleles)
    """
    print("no variant", sample['no variant'])
    print("in frame", in_frame_porc)
    print("out frame", out_frame_porc)
    print("alleles in", str(len(alleles)))
    print("alleles out", str(len(frame_shift_alleles)))
    print(out_frame_porc + in_frame_porc)
    """
    return in_frame_porc*100
