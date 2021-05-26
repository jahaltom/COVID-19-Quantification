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



#####Read in run_accession ids######
with open ("RAids.txt") as f:
        ra_ids=f.read().splitlines()

#create salmon object. If index is not present it will be created
salmon=quant.Salmon(index="human_data/salmon_index",transcriptome=Tr,threads=15)

rule all:
	input:
		expand("{wd}/results_TPM_tx.tsv",wd=DIR)
rule quant:
	output:
		quant_file="{wd}/{sample}/salmon_out/quant.sf"
	run:
		#Path to quant file
		outfile=str(output.quant_file)
		#get srrid{sample}
		srrid=outfile.split("/")[1]
		#Run Salmon on sra object and delete fastq when finished. 
		sra.SRA(srrid,directory=DIR).quant(salmon).delete_fastq()

rule merge:
	input:
		["{wd}/{sample}/salmon_out/quant.sf".format(wd=DIR,sample=s) for s in ra_ids]
	output:
		outfile="{wd}/results_TPM_tx.tsv"
	run:
		#read in file
		with open(input[0]) as f:
			thisdata=f.read().splitlines()
		thisdata.pop(0)		
        	#Run accession IDs
		names=[]
        	##Transcript IDs
		txids=[]       
		for l in thisdata:
            	#Get 1st column (TranscriptID)
			thistx=l.split('\t')[0]                    
			txids.append(thistx)
		##For TPM	
		dftpm=pd.DataFrame({'TranscriptID':txids})
        	##Make a copy for the counts
        	dfcount=pd.DataFrame({'TranscriptID':txids})
		#read files in memory
		for qf in input:
            		##current filename(RAid)
			name=qf.split('/')[1]
            		##Will become list of all RAids
			names.append(name)
            		##Get TPM
			thisdata=pd.read_csv(qf,sep='\t',usecols=[3],skiprows=0)
            		##Get Counts
            		counts=pd.read_csv(qf,sep='\t',usecols=[4],skiprows=0)
            		##Add TPM and counts to respective df with current Run accession ID as column name
			dftpm[name]=thisdata['TPM']
            		dfcount[name]=counts['NumReads']

        	#transcript TPMs. 
		df_tx=dftpm[['TranscriptID']+names].copy()
		#write to file
		df_tx.to_csv(output.outfile,sep='\t',index=False)
        
              
        	
		###Collapse transcript to respective gene names sum counts and TPMs
		
		#Gather counts and TPMs for each transcript
		df_gene_tpm=dftpm[['TranscriptID']+names].copy()	
        	df_gene_count=dfcount[['TranscriptID']+names].copy()
		
        
        	#Read in gene metadata
        	md=pd.read_csv('meta_data.tsv',sep='\t',skiprows=0) 
        	md.rename(columns={ md.columns[0]: "TranscriptID" }, inplace = True)
        	
		#Make df with transcript id and corresponding gene id
		md2=md[['TranscriptID','Gene_stable_ID']]
		##Merge TPM and count data by TranscriptID
		df_gene_tpm=md2.merge(df_gene_tpm, on=['TranscriptID'], how='right')       	      
        	df_gene_count=md2.merge(df_gene_count, on=['TranscriptID'], how='right')
		
		#Collapse so that each gene id is listed once. sum up corresponding transcript TPM and counts.
		df_gene_tpm = df_gene_tpm.groupby(['Gene_stable_ID'],as_index = False).sum()
        	df_gene_count = df_gene_count.groupby(['Gene_stable_ID'],as_index = False).sum()
		
		#Format metadata to have each gene ID once
        	del md['TranscriptID']
       		md=md.drop_duplicates()

     		
		##Merge metadata to counts and TPM
        	df_gene_tpm=md.merge(df_gene_tpm, on=['Gene_stable_ID'], how='right')        	    	
        	df_gene_count=md.merge(df_gene_count, on=['Gene_stable_ID'], how='right')
    
		#reorder so gene name is first.	
        	df_gene_tpm = df_gene_tpm[ ['Gene_name'] + [ col for col in df_gene_tpm.columns if col != 'Gene_name' ] ]
        	df_gene_count = df_gene_count[ ['Gene_name'] + [ col for col in df_gene_count.columns if col != 'Gene_name' ] ]
		
		
		
		df_gene_tpm.to_csv(DIR+'/results_TPM_gene.tsv',sep='\t',index=False)        
        	df_gene_count.to_csv(DIR+'/results_Count_gene.tsv',sep='\t',index=False)
