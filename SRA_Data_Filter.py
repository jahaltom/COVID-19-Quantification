import pandas as pd
from pandas import DataFrame
import os
from os import path


#Import TPM file.
dftpm = pd.read_csv("results_TPM_gene.tsv",sep='\t')
#Import count file.
dfcount = pd.read_csv("results_Count_gene.tsv",sep='\t')

#Import list of study_accession and run_accession IDs.
ra_list = pd.read_csv("List",sep='\t')
##Create dict that has study_accession ID as key.
d = ra_list.groupby('study_accession')['run_accession'].apply(list).to_dict()

##Will be used to get all genes that pass TPM filter. 
genes=[]

for key in d:
    #tpm will be df of all run_accession IDs and their TPM. Each Key is a single study_accession ID.
    tpm = dftpm[d[key]]
    #count will be df of all run_accession IDs and their count. Each Key is a single study_accession ID.
    count = dfcount[d[key]]
    #Metadata from TPM file. Same as count.
    metadata = dftpm[dftpm.columns[0:6]]
    ##Merge metadata with tpm and count
    tpm = pd.concat([metadata, tpm], axis=1)
    count = pd.concat([metadata, count], axis=1)
    if path.exists(key) is not True:
        os.mkdir(key)
    
    #create a median column
    tpm['median']=tpm.median(axis=1)
    #Remove genes where TPM median < 1, except for SARs genes.
    indexNames = tpm[ (tpm['median'] < 1) & (tpm['chr'] != 'SARSCOV2_ASM985889v3') ].index
    tpm.drop(indexNames , inplace=True)

    #calculate medians of median tpm dist for protein_coding, lncRNA, and EB genes
    med_med_pc=tpm.loc[tpm['Gene_type'] == 'protein_coding']['median'].median()
    med_med_lnc=tpm.loc[tpm['Gene_type'] == 'lncRNA']['median'].median()
    med_med_eb=tpm.loc[tpm['Gene_type'] == 'EB_novel']['median'].median()
    #Remove EB genes where median TPM is less than that of the med_med_lnc
    indexNames = tpm[(tpm['Gene_type'] == 'EB_novel') & (tpm['median'] < med_med_lnc) ].index
    tpm.drop(indexNames , inplace=True)
    ##Append genes that made it through filter to list
    genes.append(tpm['Gene_stable_ID'])
    #write to respective study
    tpm.to_csv(key +'/'+ key + ".TPM.tsv",sep='\t',index=False)
    count.to_csv(key +'/'+ key + ".Count.tsv",sep='\t',index=False)
    
    
##Use genes that pass TPM filter(genes) to filter file that contains TPM for all studies. 
genes = pd.concat(genes,ignore_index=True)
genes = genes.drop_duplicates()
genes = DataFrame(genes)
genes.columns=['Gene_stable_ID']
dftpm = pd.merge(dftpm,genes,on=['Gene_stable_ID'])
dftpm.to_csv("results_TPM_gene.filtered.tsv",mode='w', header=True,index=False)




