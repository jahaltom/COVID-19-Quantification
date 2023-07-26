import pandas as pd
from pandas import DataFrame



#Import TPM file.
dftpm = pd.read_csv("results_TPM_tx.tsv",sep='\t')
#Import count file.
dfcount = pd.read_csv("results_count_tx.tsv",sep='\t')


#Import metadata with samlpe IDs and corrisponding condition included.
md = pd.read_csv("metadata.tsv",sep='\t')
##Create dict that has Condition as key.
d = md.groupby('Condition')['SampleID'].apply(list).to_dict()

##Will be used to get all tx and genes that pass TPM filter.
tx=[]
gene=[]
##Median info
infoOut=[]

for key in d:
    #tpm will be df of all run_accession IDs and their TPM. Each Key is a Condition.
    tpm = dftpm[d[key]]
    #Metadata from TPM file.
    metadata = dftpm[dftpm.columns[0:37]]
    #create a median column.
    tpm['median']=tpm.median(axis=1)
    ##Merge metadata with tpm.
    tpm = pd.concat([metadata, tpm], axis=1)
    #Remove tx where TPM median < 1. Except in SARS-COV-2.
    indexNames = tpm[ (tpm['median'] < 1) & (tpm['chr'] != 'SARSCOV2_ASM985889v3')].index
    tpm.drop(indexNames , inplace=True)
    #calculate the median of medians tpm for protein_coding, lncRNA, and EB tx
    med_med_pc=tpm.loc[tpm['Gene_type'] == 'protein_coding']['median'].quantile(0.75)
    med_med_lnc=tpm.loc[tpm['Gene_type'] == 'lncRNA']['median'].quantile(0.75)
    med_med_eb=tpm.loc[tpm['Gene_type'] == 'EB_novel']['median'].quantile(0.75)
    #Remove EB tx where median TPM is less than that of the med_med_pc
    #indexNames = tpm[(tpm['Gene_type'] == 'EB_novel') & (tpm['median'] < med_med_pc) & (tpm['is54K_EB'] == False)].index
    indexNames = tpm[(tpm['Gene_type'] == 'EB_novel') & (tpm['median'] < med_med_pc)].index
    tpm.drop(indexNames , inplace=True)
    ##Append tx that made it through filter to list and record
    tx.append(tpm['TranscriptID'])
    gene.append(tpm['Gene_ID_ver'])
    tpm.to_csv(key+".TPM_EB_Filtered.tx.tsv",mode='w', sep='\t',index=False)
    info=[key,len(d[key]),"Protein Coding:",med_med_pc,"lincRNA:",med_med_lnc,"Evidence based:",med_med_eb,"Total:",len(tpm[tpm['Gene_type'] == 'EB_novel'])]
    infoOut.append(info)


##Use genes that passed the TPM filter(genes) to pull from a file that contains TPM and Counts for all conditions. Creates median info file.
tx = pd.concat(tx,ignore_index=True)
tx = tx.drop_duplicates()
tx = DataFrame(tx)
#Make it so that only the samples from metadata.tsv are included in the filtered quant output.
dftpm=pd.concat([metadata, dftpm[md["SampleID"].tolist()]], axis=1)
dfcount=pd.concat([metadata, dfcount[md["SampleID"].tolist()]], axis=1)

dftpm = pd.merge(dftpm,tx,on=['TranscriptID'])
dftpm.to_csv("results_TPM_tx.EBfiltered.tsv",mode='w', sep='\t',index=False)
dfcount = pd.merge(dfcount,tx,on=['TranscriptID'])
dfcount.to_csv("results_Count_tx.filtered.tsv",mode='w', sep='\t',index=False)
infoOut=pd.DataFrame(infoOut)
infoOut.to_csv("medianTPMInfo.tsv",mode='w', sep='\t',header=None,index=False)


#Read in gene stuff
##Use genes that passed the TPM filter(genes) to pull from a file that contains TPM and Counts for all conditions. Creates median info file.
gene = pd.concat(gene,ignore_index=True)
gene = gene.drop_duplicates()
gene = DataFrame(gene)

#Import TPM file.
dftpmG = pd.read_csv("results_TPM_gene.tsv",sep='\t')
#Import count file.
dfcountG = pd.read_csv("results_Count_gene.tsv",sep='\t')

metadata = dftpmG[dftpmG.columns[0:30]]

#Make it so that only the samples from metadata.tsv are included in the filtered quant output.
dftpmG=pd.concat([metadata, dftpmG[md["SampleID"].tolist()]], axis=1)
dfcountG=pd.concat([metadata, dfcountG[md["SampleID"].tolist()]], axis=1)

dftpmG = pd.merge(dftpmG,gene,on=['Gene_ID_ver'])
dftpmG.to_csv("results_TPM_gene.EBfiltered.tsv",mode='w', sep='\t',index=False)
dfcountG = pd.merge(dfcountG,gene,on=['Gene_ID_ver'])
dfcountG.to_csv("results_Count_gene.EBfiltered.tsv",mode='w', sep='\t',index=False)




