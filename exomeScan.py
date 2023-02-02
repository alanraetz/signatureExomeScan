# -*- coding: utf-8 -*-
"""
Extract COSMIC mutational signature (relative frequency of each nucleotide triplet - 96 entries/rows)
  - each row represents the frequency at that nucleotide triplet 
"""

# import re
import sys
import csv
import pandas as pd
import itertools
import numpy as np
import codon_table # call codonLookup(codon)
from itertools import chain
from tqdm import tqdm

debug = 4
aa_change = {}

def get_wt_codon(string):
    return string[0] + string[2] + string[6]

def get_mut_codon(string):
    return string[0] + string[4] + string[6]

#  [(i, x[i:i+2]) for i in findall('na', x)]

def findall(p, s):
    '''Yields all the positions of
    the pattern p in the string s.'''
    i = s.find(p)
    while i != -1:
        yield i      # returns found position, iterates on static function 
        i = s.find(p, i+1)  # includes overlapping matches 

def get_signature (sigType):
    
    sigType = 'SBS7a'
    sigs = pd.read_csv("COSMIC_v3.3_SBS_GRCh38.txt", delim_whitespace=True)
    freq            = sigs[["Type",sigType]]
    two_percent_min = freq[sigType].multiply(100) 
    test_min        =  two_percent_min >= 2   # returns boolean True or False
    # one_percent_min = one_percent_min if one_percent_min >= 1 else 0
    freq.insert(2,'percent',two_percent_min)
    freq.insert(3,'use_freq', test_min)

    result = freq[ freq['use_freq'] == True ]
    sigs = {} # dictionary
    for index,row in result.iterrows():
        sigs[ get_wt_codon(row.Type) ] = [get_mut_codon(row.Type), row.percent ]
        if debug > 2:
            print( str(index) + ' : ' + str(row.Type) + ' ' + str(row.percent))
            
    return sigs


def get_one_peptide (pro,pro_len,pos,aa1,aa2,peptide_radius):
    
    if pos+2 >= len(pro):
        mutant_protein = pro[0:pos] + aa1 + aa2
    elif aa2.isalpha() or aa2 == '*': # ENSG00000251608 special case internal stop codons 
        mutant_protein = pro[0:pos] + aa1 + aa2 + pro[pos+2:]
    else:
        print("protein length "+str(pro_len)+" index out of range at " + str(pos))
        
    #mut_len = len(mutant_protein)
    
    start = int(pos - peptide_radius) 
    if start < 0 :
        start = 0
    end = pos + peptide_radius
    if end > pro_len :
        end = pro_len
              
    fragment = mutant_protein[start:end]
    # frag_len = len(fragment)
    wt = pro[start:end]
        
    if debug > 2: 
        
        wt_aa = pro[pos:pos+2]
        wt1 = pro[pos:pos+1]
        wt2 = pro[pos+1:pos+2]
        print(' pos ' + str(pos) + ': ' + wt_aa + ' > ' + aa1 + aa2 + ' : ' + wt + ' => ' + fragment)
        if aa1 and wt1 and wt1 != aa1:
            aa_key = wt1 + aa1
            if aa_key in aa_change:
                aa_change[aa_key] += 1
            else: 
                aa_change[aa_key] = 1
        if aa2 and wt2 and wt2 != aa2:
            aa_key = wt2 + aa2
            if aa_key in aa_change:
                aa_change[aa_key] += 1
            else: 
                aa_change[aa_key] = 1

    return [wt,fragment] # return peptide fragment centered on mutation +/- offet
 
def get_mutant_peptides(pro,pro_len,pos,aa1,aa2,peptide_radius):

    # pos, aa1, and aa2 are lists of amino acids modified in this protein
    # loop through and generate all the modified peptide with flanking region
    mutant_pep = []
    wt_pep     = []
    all_pos    = []
    
    for position,p1,p2 in zip(pos,aa1,aa2):
        
        wt,mut = get_one_peptide(pro,pro_len,position,p1,p2,peptide_radius)
        
        all_pos.append(position)
        wt_pep.append(wt)    
        mutant_pep.append(mut)
         
    return all_pos,wt_pep,mutant_pep

       
def get_one_AA(id,nt_pos,mutation,dna,pro,pLen):
    
    offset = nt_pos % 3 # offset=0 is pos1, offset=1 is pos2, offset=2 is pos3 
    
    # handle zero as nt_pos ?
    nt_start  = nt_pos - offset
    aa_offset = int(nt_pos / 3) # protein positions start at zero!
  #  if debug > 2 :
  #      print("pos " + str(nt_pos) + ": offset = " + str(offset))
        
    if ( offset == 0 ): # special case, only one amino acid changed (at most)
        mut_aa1 = codon_table.codonLookup(mutation)
        index = aa_offset+1
        if 0 <= index < len(pro):
            mut_aa2 = pro[aa_offset+1]
        else:
            mut_aa2 = '' # special case: aa1 is the last AA in protein! 
    else:
        # test two amino acid positions +0 and +1
        # wt_aa1 = pro[aa_offset]
        # wt_aa2 = pro[aa_offset+1]
        # get wt codon sequence at the mutation position
        wt_codon1 = dna[nt_start:nt_start+3]
        wt_codon2 = dna[nt_start+3:nt_start+6]    
        codon1 = wt_codon1
        codon2 = wt_codon2
        
        # map mutation into two amino acid reading frame slots
        if ( offset == 1 ):
            codon1 = codon1[0:1] + mutation[0:2]
            codon2 = mutation[2:3] + codon2[1:3]  
        elif ( offset == 2 ):
            codon1 = codon1[0:2] + mutation[0:1]
            codon2 = mutation[1:3] + codon2[2:3]  
                
        # Ensembl transcript is longer than protein
        #mut_aa1 = codon_table.codonLookup(codon1)
        if len(codon1) == 3:# default case
            mut_aa1 = codon_table.codonLookup(codon1) 
        else:  # ENSG00000198727 special case? 
            mut_aa1 = '_' # stop codon
          
        # Ensembl transcript is longer than protein
        if len(codon2) == 3:# default case
            mut_aa2 = codon_table.codonLookup(codon2) 
        else:  # ENSG00000198727 special case? 
            mut_aa2 = '_' # stop codon
        
    if debug > 4 :
        print("AA " + str(aa_offset) + " = " + mut_aa1 + mut_aa2)
        
    if ( mut_aa1 == '_' or mut_aa2 == '_' ):
        return(0,0,0) # stop codon, so skip this mutation
        
    elif ( mut_aa1 == 'N' or mut_aa2 == 'N' ): 
    
        print('SKIPPING '+str(id)+' nt position '+str(nt_pos)+' = '+mutation)
        return(0,0,0) # codonTable.codonLookup failed (skip this mutation)
        
    elif ( mut_aa1 == pro[aa_offset] and mut_aa2 == pro[aa_offset+1] ):
        return(0,0,0) # silent mutations
    else:
        return (aa_offset,mut_aa1,mut_aa2)
            
def get_mutant_AAs(id,nt_pos,mutation,dna,pro,pLen):
    
    ''' 
    translate a single instance of a mutation into the resulting amino acid changes (0,1,or 2) 
    '''
    # sanity check: DNA minus stop codon / 3 should equal protein length
    dna_codons = float((len(dna))/3)
    if( dna_codons != pLen ):
        print(id + " dna = " + str(dna_codons) + ", does not equal protein = " + str(pLen))
        # sys.exit(1)
     
    aa_offset = []
    mut_aa1   = []
    mut_aa2   = []  
    for pos in nt_pos: # this is the list of positions where the 3 nt 'mutation' occurs
        
        (off,aa1,aa2) = get_one_AA(id,pos,mutation,dna,pro,pLen)
        if ( off ):
            aa_offset.append(off)
            mut_aa1.append(aa1)
            mut_aa2.append(aa2)
            
    return(aa_offset,mut_aa1,mut_aa2)


def signature_search(UV_sigs,dna):

    # gene_size = length(dna)
    mutation_list = []
    start_positions = []
    freqs = []
    for wt_seq,values in UV_sigs.items():    # for key, value in a_dict.items():    
        
        # string search through gene sequence dna to get all instances of wt_seq
        # note: matches returned without respect to amino acid reading frame 
        matches = [ i for i in findall(wt_seq,dna)] 
        if len(matches) > 0:
            start_positions.append(matches)
            mutation_list.append( values[0] ) # mutated DNA sequence
            freqs.append( values[1] ) # frequency this sequence appears

    return start_positions,mutation_list,freqs

def get_signature_peptides (signature,mutation_signature):
    
    genome_size = len(all_genes)
     
    for id,dna in tqdm(all_genes.items(), total = genome_size):
    # for id,dna in all_genes.items():

        #if id != 'ENSG00000251608': # internal stop codons
        #    continue
        if id == "ENSG00000283496": # skip, bad # "ENSG00000198727":
            continue
        (nt_positions,mutations,frequencies) = signature_search(mutation_signature,dna)
        pro = all_proteins[id]
        pLen = len(pro)
        gene_mut_count = 0
        gene_mutant_peptides = {}
        pepfile = open('exomeScanPeptides.fasta','a') # append (re-save file each gene)

        for mut_pos,mut_codon,freq in zip(nt_positions,mutations,frequencies):
            
            if debug > 2: 
                print('found signature '+ mut_codon + ' at pos ' + str(mut_pos))
                
            (aa_pos,aa1,aa2) = get_mutant_AAs(id,mut_pos,mut_codon,dna,pro,pLen)
            
            if (aa_pos): # has at least one entry              18 = max +/- offset from mutation
                    
                pos,wt,mut = get_mutant_peptides(pro,pLen,aa_pos,aa1,aa2,18)
    
                for position,wt_pep,mut_pep in zip(pos,wt,mut):

                    gene_mutant_peptides[id +'_'+ str(position) + '_wt']  = wt_pep
                    gene_mutant_peptides[id +'_'+ str(position) + '_mut'] = mut_pep
                    gene_mut_count += 1
                    if debug > 2:
                        print(id+" pos "+str(position)+" : "+ wt_pep +' > '+ mut_pep )                                                
             
        if debug > 1:
            
            print(id+" has "+str(gene_mut_count)+" "+signature+" signature mutations")
            for mut,mutCount in aa_change.items(): # global variable
                 print('mutation '+mut+' : '+str(mutCount)+' instances')
                  
        for name,peptide_seq in gene_mutant_peptides.items():
            pepfile.write(">" + name +"\n")  # fasta sequence name line
            pepfile.write(peptide_seq + "\n") # fasta sequence
      
        pepfile.close()
        


    # TODO: return signature frequency for each peptide
                
    return ()    

# _________START OF MAIN PROGRAM______________

signature = 'SBS7a'    
            
UV_sig = get_signature(signature) # must be header name in COSMIC SBS file! 

all_genes    = pd.read_csv('human-genes-dna.csv',header=None, index_col=0, squeeze=True).to_dict()
all_proteins = pd.read_csv('human-proteins.csv', header=None, index_col=0, squeeze=True).to_dict()

get_signature_peptides(signature,UV_sig)
    
#with open('signature_peptides_' + signature + '.txt', 'w') as csvFile:
#    csvFile.DictWriter.writerows([signature_peptides,sig_frequency])


#[0, 5, 10, 15]
