import pandas as pd
from pandas import DataFrame
import os
from os import path


#Import TPM file.
dftpm = pd.read_csv("results_TPM_gene.tsv",sep='\t')
#Import count file.
dfcount = pd.read_csv("results_Count_gene.tsv",sep='\t')


#Import list of study_accession and run_accession IDs.
ra_list = pd.read_csv("Study.tsv",sep='\t')
##Create dict that has study_accession ID as key.
d = ra_list.groupby('study_accession')['run_accession'].apply(list).to_dict()

for key in d:
    #tpm will be df of all run_accession IDs and their TPM. Each Key is a single study_accession ID.
    tpm = dftpm[d[key]]
    #count will be df of all run_accession IDs and their count. Each Key is a single study_accession ID.
    count = dfcount[d[key]]
    #Metadata from TPM file. Same as count.
    metadata = dftpm[dftpm.columns[0:29]]
    

    
    ##Merge metadata with tpm and count
    tpm = pd.concat([metadata, tpm], axis=1)
    count = pd.concat([metadata, count], axis=1)
    if path.exists(key) is not True:
        os.mkdir(key)


    #write to respective study
    tpm.to_csv(key +'/'+ key + ".TPM.tsv",sep='\t',index=False)
    count.to_csv(key +'/'+ key + ".Count.tsv",sep='\t',index=False)
