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
import codon_table # call codonLookup(codon)

def get_wt_codon(string):
    return string[0] + string[2] + string[6]

def get_mut_codon(string):
    return string[0] + string[4] + string[6]

def shift1(List):
    value = List[0]
    new_list = List[1: ]
    return (value,new_list)

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
            
    return sigs


def get_mutant_peptide (pro,pro_len,pos,aa1,aa2,peptide_radius):
    
    mutant_protein = pro[0:pos] + aa1 + aa2 + pro[pos+2:]
    
    start = int(pos - peptide_radius) 
    if start < 0 :
        start = 0
    end = pos + peptide_radius
    if end > pro_len :
        end = pro_len

    return mutant_protein[start:end] # return peptide fragment centered on mutation +/- offet
        
            
def get_mutant_AAs(id,nt_pos,mutation,dna,pro,pLen):
    
    ''' 
    translate a single instance of a mutation into the resulting amino acid changes (0,1,or 2) 
    '''
    # sanity check: DNA minus stop codon / 3 should equal protein length
    dna_codons = float((len(dna-3))/3)
    if( dna_codons != pLen ):
        sys.exit(id + " dna = " + dna_codons + ", does not equal protein = " + pLen)
        
    offset = nt_pos % 3 # offset=0 is pos1, offset=1 is pos2, offset=2 is pos3 
    
    # handle zero as nt_pos ?
    nt_start  = nt_pos - offset
    aa_offset = int(nt_pos / 3) # protein positions start at zero!
    
    if ( offset == 0 ): # special case, only one amino acid changed (at most)
        mut_aa1 = codon_table.codonLookup[mutation]
        mut_aa2 = pro[aa_offset+1]
    else:
        # test two amino acid positions +0 and +1
        # wt_aa1 = pro[aa_offset]
        # wt_aa2 = pro[aa_offset+1]
        # get wt codon sequence at the mutation position
        wt_codon1 = dna[nt_start:nt_start+3]
        wt_codon2 = dna[nt_start+3:nt_start+7]
        
        codon1 = wt_codon1
        codon2 = wt_codon2
        
        # map mutation into two amino acid reading frame slots
        if ( offset == 1 ):
            codon1 = codon1[0:1] + mutation[0:2]
            codon2 = mutation[2:3] + codon2[1:3]  
        elif ( offset == 2 ):
            codon1 = codon1[0:2] + mutation[0:1]
            codon2 = mutation[1:3] + codon2[2:3]  
        else:
            sys.exit(id + ": bad offset = " + offset + " for nt pos " + nt_start)
            
        mut_aa1 = codon_table.codonLookup(codon1)
        mut_aa2 = codon_table.codonLookup(codon2)
        
    if ( mut_aa1 == '_' or mut_aa2 == '_' ):
        return() # stop codon, so skip this mutation
    elif ( mut_aa1 == pro[aa_offset] and mut_aa2 == pro[aa_offset+1] ):
        return() # silent mutations
    else:
        return (aa_offset,mut_aa1,mut_aa2)

def signature_search(UV_sigs,dna):

    # gene_size = length(dna)
    mutation_list = []
    start_positions = []
    freqs = []
    for wt_seq,values in UV_sigs.items:        
        
        # string search through gene sequence dna to get all instances of wt_seq
        # note: matches returned without respect to amino acid reading frame 
        if start_positions.append( [ i for i in findall(wt_seq,dna)] ):
            mutation_list.append( values[0] ) # mutated DNA sequence
            freqs.append( values[1] ) # frequency this sequence appears

    return start_positions,mutation_list,freqs

def get_signature_peptides (mutation_signature):
    
    mutant_peptide_list = []
    mutant_frequencies = []
    
    for id,dna in all_genes.items():

        (nt_positions,mutations,frequencies) = signature_search(mutation_signature,dna)
        pro = all_proteins[id]
        pLen = len(pro)

        for mut_pos in nt_positions:
            
            (mutation,mutations) = shift1(mutations)
            (aa_pos,aa1,aa2) = get_mutant_AAs(id,mut_pos,mutation,dna,pro,pLen)
            
            if (aa_pos):        # 15 = max +/- offset from mutation
                mutant_peptide_list.append(get_mutant_peptide(pro,pLen,aa_pos,aa1,aa2,15))
                mutant_frequencies.append(frequencies)
            else:
                next()
                
    # TODO: return signature frequency for each peptide
                
    return (mutant_peptide_list,mutant_frequencies)    

# _________START OF MAIN PROGRAM______________

signature = 'SBS7a'    
            
UV_sig = get_signature(signature) # must be header name in COSMIC SBS file! 

all_genes    = pd.read_csv('human-genes-dna.csv', header=None, index_col=0, squeeze = True).to_dict()
all_proteins = pd.read_csv('human-proteins.csv',  header=None, index_col=0, squeeze = True).to_dict()

(signature_peptides,sig_frequency) = get_signature_peptides(UV_sig)
    
with open('signature_peptides_' + signature + '.txt', 'w') as csvFile:
    csvFile.DictWriter.writerows([signature_peptides,sig_frequency])


#[0, 5, 10, 15]
