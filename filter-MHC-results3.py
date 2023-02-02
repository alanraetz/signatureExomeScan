# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:49:47 2022

@author: alanr
"""
wt = {}
mut = {}
full_row = {}
import csv
with open('mhcflurry-results.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        idx,name,pos,peptide,n_flank,c_flank,sample,affinity,best_allele,affpct,pro_score,pre_score,prepct = row
        full_row[name] = row # save the row for the next loop
        
        ensg_pos = name.find('_wt')
        if ( ensg_pos != -1 ):
            wt[name[0:ensg_pos]] = affinity  # create a WT dict entry with the ENSG name (without suffix) is the key
        else:
            ensg_pos = name.find('_mut')
            if ( ensg_pos != -1 ):
                mut[name[0:ensg_pos]] = affinity # create a MUT dict entry with the ENSG name (without suffix) is the key

f = open('filtered-mhc-results_1008.csv', 'w')
writer = csv.writer(f)
count = 0
for ensg in mut: # we only care about low affinity mutants

    if ensg in wt:
        
        if mut[ensg] < wt[ensg]:
            diff = float(wt[ensg]) - float(mut[ensg])
            count +=1
            print(ensg + " wt > mut :" + str(wt[ensg])+'_> '+str(mut[ensg]))
            row2 = full_row[ensg+"_mut"]
            affinity = row2[7]
            percent = (diff / (diff+float(affinity)))*100
            row2.append(int(diff))
            row2.append(int(percent))
            writer.writerow( row2 ) # full_row[ensg+"_mut"] ) 
            count += 1

f.close()

print("Found "+str(count)+" mutant peptides with higher affinity")
            