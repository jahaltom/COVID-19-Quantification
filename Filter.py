import pandas as pd
from pandas import DataFrame
import os
from os import path


#Import TPM file.
dftpm = pd.read_csv("results_TPM_gene.tsv",sep='\t')
#Import count file.
dfcount = pd.read_csv("results_Count_gene.tsv",sep='\t')

#Import a list of samlpe IDs and corrisponding tissue. 
ra_list = pd.read_csv("List.tsv",sep='\t')
##Create dict that has tissue as key.
d = ra_list.groupby('Tissue')['SampleID'].apply(list).to_dict()

##Will be used to get all genes that pass TPM filter. 
genes=[]
##Median info
med_info=[]

for key in d:
    #tpm will be df of all run_accession IDs and their TPM. Each Key is a tissue.
    tpm = dftpm[d[key]]
    #count will be df of all run_accession IDs and their count. Each Key is a tissue. 
    count = dfcount[d[key]]
    #Metadata from TPM file. Same as count.
    metadata = dftpm[dftpm.columns[0:20]]
    
    #create a median column
    tpm['median']=tpm.median(axis=1)
    
    ##Merge metadata with tpm and count
    tpm = pd.concat([metadata, tpm], axis=1)
    count = pd.concat([metadata, count], axis=1)
    
    #Remove genes where TPM median < 1
    indexNames = tpm[ (tpm['median'] < 1)].index
    tpm.drop(indexNames , inplace=True)
    
    #calculate the median of medians tpm for protein_coding, lncRNA, and EB genes
    med_med_pc=tpm.loc[tpm['Gene_type'] == 'protein_coding']['median'].median()
    med_med_lnc=tpm.loc[tpm['Gene_type'] == 'lncRNA']['median'].median()
    med_med_eb=tpm.loc[tpm['Gene_type'] == 'EB_novel']['median'].median()
    info=[key,"Protein Coding:",med_med_pc,"lincRNA:",med_med_lnc,"Evidence based:",med_med_eb]
    med_info.append(info)
    
    #Remove EB genes where median TPM is less than that of the med_med_pc
    indexNames = tpm[(tpm['Gene_type'] == 'EB_novel') & (tpm['median'] < med_med_pc) & (tpm['is54K_EB'] == False)].index
    tpm.drop(indexNames , inplace=True)
    ##Append genes that made it through filter to list
    genes.append(tpm['Gene_stable_ID'])
    
    
       
##Use genes that passed the TPM filter(genes) to pull from a file that contains TPM and Counts for all tissues 
genes = pd.concat(genes,ignore_index=True)
genes = genes.drop_duplicates()
genes = DataFrame(genes)
genes.columns=['Gene_stable_ID']
dftpm = pd.merge(dftpm,genes,on=['Gene_stable_ID'])
dftpm.to_csv("results_TPM_gene.filtered.tsv",mode='w', sep='\t',header=True,index=False)
dfcount = pd.merge(dfcount,genes,on=['Gene_stable_ID'])
dfcount.to_csv("results_Count_gene.filtered.tsv",mode='w', sep='\t',header=True,index=False)

med_info=pd.DataFrame(med_info)
med_info.to_csv("medianInfo.tsv",mode='w', sep='\t',header=None,index=False)
