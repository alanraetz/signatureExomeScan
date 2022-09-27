#! /usr/bin/env/ python 
""" 
# fasta to tsv converter 
# usage: fasta2csv.py input.fasta > output.csv 
# output: dump to standard stdout 
# Dev: __author__ = 'aung' 
# Date: 2013 03 16 
""" 
import sys 
from Bio import SeqIO 

def fasta2csv(fasta_file,output_file):
    # fasta_file = sys.argv[1] 
    with open(output_file,"w") as out:
        for seq_record in SeqIO.parse(fasta_file, "fasta"): 
            out.write(str(seq_record.id)+','+str(seq_record.seq)+'\n') 

output_file = 'test.csv'
fasta_file = sys.argv[1]
fasta2csv(fasta_file,output_file)


