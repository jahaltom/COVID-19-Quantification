import pandas as pd



#Import TPM file.
df_gene_tpm = pd.read_csv("results_TPM_gene.tsv",sep='\t')


#####Filter TPM
#create a median column
df_gene_tpm['median']=df_gene_tpm.median(axis=1)
#Remove genes where TPM median < 1, except for SARs genes. 
indexNames = df_gene_tpm[ (df_gene_tpm['median'] < 1) & (df_gene_tpm['chr'] != 'SARSCOV2_ASM985889v3') ].index
df_gene_tpm.drop(indexNames , inplace=True)

#calculate medians of median tpm dist for protein_coding, lncRNA, and EB genes
med_med_pc=df_gene_tpm.loc[df_gene_tpm['Gene_type'] == 'protein_coding']['median'].median()
med_med_lnc=df_gene_tpm.loc[df_gene_tpm['Gene_type'] == 'lncRNA']['median'].median()
med_med_eb=df_gene_tpm.loc[df_gene_tpm['Gene_type'] == 'EB_novel']['median'].median()
#Remove EB genes where median TPM is less than that of the med_med_lnc
indexNames = df_gene_tpm[(df_gene_tpm['Gene_type'] == 'EB_novel') & (df_gene_tpm['median'] < med_med_lnc) ].index
df_gene_tpm.drop(indexNames , inplace=True)

df_gene_tpm.to_csv('results_TPM_gene.tsv',sep='\t',index=False)        
