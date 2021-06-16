import pandas as pd


#####################Gene Level

#Read in metadata from phylostrata. Each gene is represented by 1 line, transcripts are merged. Read in general metadata. 
strata_md=pd.read_csv("gene_level_metadata.tsv",sep='\t')
md=pd.read_csv("meta_data.tsv",sep='\t')

#Extract SARS transcripts and remove transcript column.
sarsmd=md[md['Gene_stable_ID'].str.contains('ENSSAST')]
#for later
sarsmd2=sarsmd.copy()
sarsmd.drop('tid', axis=1, inplace=True)

#Remove SARS transcripts. Only Gencode 36 and EB genes have phylostrata. 
md=md[~md['Gene_stable_ID'].str.contains('ENSSAST')]
#Unwanted columns
md = md.drop(['Gene_name', 'strand','Gene_description','tid'], axis=1)
#Change EB Gene_stable_ID to match the current version. 
md['Gene_stable_ID'] = md['Gene_stable_ID'].str.replace("mikado", "EB")
md=md.drop_duplicates()

#Merge two metadatas together to add phylostrata info to other. 
md = md.merge(strata_md, on=['Gene_stable_ID','chr'])

#Merge and write
gene_md=[md,sarsmd] 
gene_md= pd.concat(gene_md)
##Once the snakefile is fun, this file will be compressed to the gene level##
gene_md.to_csv('Gene_level_metadata.tsv',sep='\t',index=False)        



#################Transcript level
##Read in protein seq's for Human EB genr, gencode gene, and Sars gene transcripts, and transcript phylostrata metadata. Rename columns. 
eb=pd.read_csv("EB",sep='\t',header=None)
eb.columns =['TranscriptID', 'seq']
human=pd.read_csv("Human",sep='\t',header=None)
human.columns =['TranscriptID', 'seq']
sars=pd.read_csv("sars",sep='\t',header=None)
sars.columns =['TranscriptID', 'seq']
sarsmd2.rename(columns={'tid': 'TranscriptID'}, inplace=True)
tx=pd.read_csv("transcript_level_metadata.tsv",sep='\t')
tx.rename(columns={'tid': 'TranscriptID'}, inplace=True)



#Get those transcripts where there was no protein coded for
nc_tx= tx[tx['final_strata'].isnull()]
#Merge metadata with EB
eb_tx = tx.merge(eb, on=['TranscriptID'])
#Merge with Gencode
human_tx = tx.merge(human, on=['TranscriptID'])
#Merge with Sars
sars_tx = sarsmd2.merge(sars, on=['TranscriptID'])

#Merge and write
tx_md=[eb_tx,human_tx,nc_tx,sars_tx]
tx_md = pd.concat(tx_md)
tx_md.to_csv("Transcript_level_metadata.tsv",mode='w', header=True,index=False,sep='\t')
