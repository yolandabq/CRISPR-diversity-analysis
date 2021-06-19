# -*- coding: utf-8 -*-
"""
Created on Wed May 19 18:01:58 2021

@author: yolib

Functions related to microhomology. 

"""


def MH_detection(seq, allele_cigar, PC = 30): 
    
    # Function to obtain the MH sequence in a deletion mediated by MH. 
    # The input is the target sequence and the cigar of the deletion we want to analyze. 
    # PC is the position of the cleavage. 

    MH = True
    MH_seq_list = []
    MH_seq = ''
    
    
    pos = int(allele_cigar.split(":")[0])
    size = int(allele_cigar.split(":")[1][:-1])
    
    if size <= 2: 
        return MH_seq
    
    pos_inicio = PC + pos - 1 
    
    while MH: # We keep looking for the MH until we don't find it (MH = False)
       
        
        try:
            if seq[pos_inicio] == seq[pos_inicio+size]: # If the starting position is the same than the start of the indel
                
                MH_seq_list.append(seq[pos_inicio])
                pos_inicio = pos_inicio - 1
                
                
            else: 
                MH = False
                
        except IndexError: # in case it goes out of the sequence
            MH = False
    
    if len(MH_seq_list) >= 2: # if the MH is >= 2, we merge the nts
        MH_seq = ''.join([str(elem) for elem in MH_seq_list])
    
    else: 
        pos_inicio = PC + pos # If we can't find MH, we verify the other direction 
        MH = True
        MH_seq_list = []
        
        while MH == True:
            try: 
                if seq[pos_inicio] == seq[pos_inicio+size]: 
                    MH_seq_list.append(seq[pos_inicio])
                    pos_inicio = pos_inicio + 1
                    
                else: 
                    
                    MH = False
                    
            except IndexError: 
                MH = False
                
        if len(MH_seq_list) >= 2:
            MH_seq = ''.join([str(elem) for elem in MH_seq_list])
        
    
    return MH_seq


def n_MH_alleles(sequence, list_deletion_cigar, PC = 30): 
    
    # Function to obtain the number of MH deletions of a mutational profile. 
    # The input is the sequence and the list of deletion cigars. 
    # PC is the position of the cleavage. 
    
    MH_list = []
    
    for deletion in list_deletion_cigar: 
        MH = MH_detection(sequence, deletion, PC = PC)
        
        if MH != '': 
            MH_list.append(MH)
            
    return len(MH_list)

def all_MH_cigar(seq, cigar_deletion_list, PC = 30):
    
    # Function that, from the list of deletions, finds the deletions mediated by MH. 
    # The output is the list of cigars produced by MH. 
    # The input is the sequence and the list of deletion cigars.
    # PC is the position of the cleavage. 
    
    MH_mut_list = []
    long_seq = len(seq)
    
    for allele in cigar_deletion_list: 
            
            # We verify if the deletion overlap with the DSB 
            
            in_PC = True
            MH_sequence = ''
            
            if int(allele.split(":")[0]) < 0 and (int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) >=0 and (PC + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= long_seq: 
                MH_sequence = MH_detection(seq, allele)
                
            #elif int(allele.split(":")[0]) == 1 and (PC + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= long_seq: # si la posición inicial es 1, también pasa por el punto de corte
                
                #MH_sequence = MH_detection(seq, allele)
                
                
            elif int(allele.split(":")[0]) > 1 and (PC + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= long_seq: 
                #if the starting position is greater than 1, we try an alternative alignment (if possible),
                # and we check again if it passes through the DSB.
                pos_inic_indel = PC + int(allele.split(":")[0]) - 2
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
                    
            elif int(allele.split(":")[0]) < 0 and (int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) < 0 and (PC + int(allele.split(":")[0]) + int(allele.split(":")[1][:-1])) <= long_seq and int(allele.split(":")[0]) >= -long_seq: 
            # if the starting position starts before the DSB but also ends before the DSB, we also try to realign
            # (but the indices change as to whether the indel is on the right)
                pos_inic_indel = PC + int(allele.split(":")[0]) 
                alt_alin = True
                alt_pos = 0
                while alt_alin:
                    if seq[pos_inic_indel] == seq[pos_inic_indel + int(allele.split(":")[1][:-1])]: 
                        alt_pos = alt_pos + 1
                        pos_inic_indel = pos_inic_indel + 1
            
                    else: 
                        alt_alin = False
                        
                alt_allele_cigar = int(allele.split(":")[0]) + alt_pos
                alt_allele_cigar_str = str(alt_allele_cigar) + ":" + allele.split(":")[1][:-1] + "D"
                
                if int(alt_allele_cigar_str.split(":")[0]) < 0 and (int(alt_allele_cigar_str.split(":")[0]) + int(alt_allele_cigar_str.split(":")[1][:-1])) >=0: 
                    MH_sequence = MH_detection(seq, alt_allele_cigar_str)
                    
                #elif int(alt_allele_cigar_str.split(":")[0]) == 1: 
                    
                #    MH_sequence = MH_detection(seq, alt_allele_cigar_str)
        
                else: 
                    
                    in_PC = False
                    
            if MH_sequence == '' and in_PC == True: 
                
                #NO_MH_mut_list.append(allele)
                pass
            
            elif MH_sequence != '' and in_PC == True: 
                MH_mut_list.append(allele)
                #MH_seq_list.append(MH_sequence)
    
    return MH_mut_list




def MH_patterns(sequence, PC = 29, file = '..\\MH\\file_MH_score.txt', save = True, output = False):

    
    # Function to obtain all the possible deletions mediated by MH from a target sequence.
    # Input: 
        # Target sequence. 
        # File: name of the file we want to create. 
        # save: Boolean to choose if we want to save the file. 
        # output: Boolean to choose if we want the dataframe as output. 
        # PC: Cleavage position. 

    # The output is a file with the next information: 
        # MH sequence, MH length, the start position of the MH at the righ of the DSB, the start position of the MH at the left of the DSB, 
        # the sequence of the deletion, the length of the deletion, the cigar of the allele (the deletion), the % GC, the % AT and if the deletion is 
        # out of frame or not (True/False). 
        
    import pandas as pd 
    
    left = PC

    right = len(sequence) - left
    
    df = pd.DataFrame(columns = ["MH_sequence", "MH_length", "Start_left", "Start_right", "Deletion_seq", "Deletion_length", "Allele", "GC", "AT", "Out_of_frame"])
    
   # with open(file, 'w') as f:
        
   #     f.write("MH_sequence\tMH_length\tStart_left\tStart_right\tGC\tAT\tDeletion_seq\tDeletion_length\tOut_frame\tAllele\n")
        
    for p_l in range(2,left + 2)[::-1]: 
        
        #print(p_l)
        for MH_len in range(2, p_l+1): 
            #print("MH ------> " + str(MH_len))
            for p_r in range (left + 1, left + right - MH_len + 1): 
                #print(p_r)
                #print("left ---> " +  sequence[p_l - MH_len : p_l])
                #print("right ---> " + sequence[p_r: p_r + MH_len])
                #print(sequence[p_l - MH_len: p_r])
                #print(len(sequence[p_l - MH_len: p_r]))
                if sequence[p_l - MH_len : p_l] == sequence[p_r : p_r + MH_len]: 
                    MH_seq = sequence[p_l - MH_len : p_l]
                    long_MH = len(MH_seq)
                    deletion_seq = sequence[p_l - MH_len: p_r]
                    long_deletion = len(deletion_seq)
                    
                    p_inic_left = p_l - long_MH - 1
                    
                    
                    # We do a realignment (we send all the MHs to the left) so that they all have the same reference,
                     # because if not, we can have some repeating MH patterns. We change the start positions of the MH and the cigar allele string.
                     
                    while sequence[p_inic_left] == sequence[p_inic_left + long_deletion] and p_inic_left >= 0:  
                       
                        MH_seq = sequence[p_inic_left] + MH_seq
                        p_inic_left = p_inic_left - 1
                        
                    allele_str = str(p_inic_left - PC) + ":" + str(len(deletion_seq)) + "D"
                    GC = MH_seq.count("G") + MH_seq.count("C")
                    long_MH = len(MH_seq)
                    
                    if len(deletion_seq) % 3 == 0: 
                        out_frame = False
                        
                    else: 
                        out_frame = True
                        
                    
                        
                    #print("left ---> " +  sequence[p_l - MH_len : p_l])
                    #print("right ---> " + sequence[p_r: p_r + MH_len])
                    #print(sequence[p_l - MH_len: p_r])
                    #print(len(sequence[p_l - MH_len: p_r]))
                    
                    df=df.append({'MH_sequence' : MH_seq, "MH_length" : long_MH, "Start_left" : p_inic_left + 2, 
                                  "Start_right" : p_inic_left + long_deletion + 2, "Deletion_seq": sequence[p_inic_left + 1: p_inic_left + long_deletion + 1],
                                  "Deletion_length" : long_deletion,
                                  "Allele" : allele_str, "GC" : GC, "AT" : long_MH - GC, "Out_of_frame" : out_frame}, ignore_index=True)
                    
                    # we add the cigar allele because sometimes the results given by the MH are the same ("duplicates")
                    
                        #f.write(MH_seq + "\t" + str(len(MH_seq)) + "\t" + str(p_l-len(MH_seq)) + "\t" + str(p_r) + "\t" + str(GC) + "\t" + str(len(MH_seq) - GC) + "\t" + deletion_seq + "\t" + str(len(deletion_seq)) + "\t" + str(out_frame) + "\t" + allele_str + "\n")

    df = df.drop_duplicates()
    
    # we save the csv
    if save: 
        df.to_csv(file, sep = '\t', index = False)
    
    if output: 
        return df
    

def pattern_score(file = '..\\MH\\file_MH_score.txt'): 
    
    # Function to compute the pattern score of each possible deletion mediated by MH. 
    # The input is the file created by the function "MH_patterns"
    # The output is a dataframe that contains the information from the input file and a new column with the pattern scores. 
    
    import pandas as pd 
    import numpy as np
    
    df = pd.read_csv(file, sep = "\t")
    
    df["Pattern_Score"] = ((df["AT"] + (df["GC"]*2))*np.exp(-(df["Deletion_length"]/20)))*100
    # Afterwards, we sort and remove the duplicates that give the same result (we keep the ones with the highest scores)
    
    df = df.sort_values(by = "Pattern_Score", ascending = False)
    df = df.drop_duplicates(['Allele'], keep='first')

    return df
    
