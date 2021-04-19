import yaml
import sys
import os
from pyrpipe import sra,quant
import pandas as pd


####Read config#####
configfile: "config.yaml"
DIR = config['DIR']
THREADS=config['THREADS']
Tr=config['transcriptome']



#####Read SRR ids######
with open ("SRRids.txt") as f:
	SRR=f.read().splitlines()

#create salmon object. If index is not present it will be created
salmon=quant.Salmon(index="human_data/salmon_index",transcriptome=Tr,threads=15)

rule all:
	input:
		expand("{wd}/results_TPM_tx.tsv",wd=DIR)
rule quant:
	output:
		quant_file="{wd}/{sample}/salmon_out/quant.sf"
	run:
		outfile=str(output.quant_file)
		srrid=outfile.split("/")[1]
		sra.SRA(srrid,directory=DIR).quant(salmon).delete_fastq()

rule merge:
	input:
		["{wd}/{sample}/salmon_out/quant.sf".format(wd=DIR,sample=s) for s in SRR]
	output:
		outfile="{wd}/results_TPM_tx.tsv"
	run:
		#read in file
		with open(input[0]) as f:
			thisdata=f.read().splitlines()
		thisdata.pop(0)		
        	#SRR
		names=[]
        	##Transcript ID
		txids=[]       
		for l in thisdata:
            	#Get 1st column (TranscriptID)
			thistx=l.split('\t')[0]                    
			txids.append(thistx)
			
		df=pd.DataFrame({'TranscriptID':txids})
        	##Make a copy
        	dfcount=pd.DataFrame({'TranscriptID':txids})
		#read files in memory
		for qf in input:
            		##current filename(SRRid)
			name=qf.split('/')[1]
            		##Will become list of all SRRids
			names.append(name)
            		##Get TPM
			thisdata=pd.read_csv(qf,sep='\t',usecols=[3],skiprows=0)
            		##Get Counts
            		counts=pd.read_csv(qf,sep='\t',usecols=[4],skiprows=0)
            		##Add TPM and counts to respective df with current SRR id as column name
			df[name]=thisdata['TPM']
            		dfcount[name]=counts['NumReads']

        	#transcript TPMs. 
		df_tx=df[['TranscriptID']+names].copy()
		#write to file
		df_tx.to_csv(output.outfile,sep='\t',index=False)
        
              
        	
		
		#Metadata TPMs
		df_gene=df[['TranscriptID']+names].copy()
	
        
        	df_gene_count=dfcount[['TranscriptID']+names].copy()
		
        
        	#add gene names
        	md=pd.read_csv('meta_data.tsv',sep='\t',skiprows=0) 
        	md.rename(columns={ md.columns[0]: "TranscriptID" }, inplace = True)
        	
		md2=md[['TranscriptID','Gene_stable_ID']]
		df_gene=md2.merge(df_gene, on=['TranscriptID'], how='right')       	
        	##Add to counts df too
        	df_gene_count=md2.merge(df_gene_count, on=['TranscriptID'], how='right')
		
		df_gene = df_gene.groupby(['Gene_stable_ID'],as_index = False).sum()
        	df_gene_count = df_gene_count.groupby(['Gene_stable_ID'],as_index = False).sum()
		
        	del md['TranscriptID']
       		md=md.drop_duplicates()

     		
		
        	df_gene=md.merge(df_gene, on=['Gene_stable_ID'], how='right')
        	
        	##Add to counts df too
        	df_gene_count=md.merge(df_gene_count, on=['Gene_stable_ID'], how='right')
        
        
        	
        
        
		#reorder	
        	df_gene = df_gene[ ['Gene_name'] + [ col for col in df_gene.columns if col != 'Gene_name' ] ]
        	df_gene_count = df_gene_count[ ['Gene_name'] + [ col for col in df_gene_count.columns if col != 'Gene_name' ] ]
		#if first col is empty
		
		df_gene.to_csv(DIR+'/results_TPM_gene.tsv',sep='\t',index=False)
        
        	df_gene_count.to_csv(DIR+'/results_Count_gene.tsv',sep='\t',index=False)
