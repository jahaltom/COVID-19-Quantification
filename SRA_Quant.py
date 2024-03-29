import yaml
import sys
import os
from pyrpipe import sra,quant,qc
import pandas as pd


#Read in config file
configfile: "config.yaml"
DIR = config['DIR']
THREADS=config['THREADS']
Tr=config['Transcriptome']


#creates a trim_galore object.
kwargs={'--fastqc' : ''}
trim_galore=qc.Trimgalore(threads=4,**kwargs)

#Read in run_accession ids
with open ("ids.txt") as f:
        ra_ids=f.read().splitlines()

#create salmon object. If index is not present it will be created
salmon=quant.Salmon(index=Tr.split("/")[0]+"/salmon_index",transcriptome=Tr,threads=THREADS)

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
		sra.SRA(srrid,directory=DIR).trim(trim_galore).quant(salmon).delete_fastq()
rule merge:
        input:
                ["{wd}/{sample}/salmon_out/quant.sf".format(wd=DIR,sample=s) for s in ra_ids]
        output:
                "{wd}/results_TPM_tx.tsv",
                "{wd}/results_count_tx.tsv"
        run:
                #read in quant file
                with open(input[0]) as f:
                        thisdata=f.read().splitlines()
                thisdata.pop(0)
                #Run accession IDs
                names=[]
                #Transcript and gene IDs
                txids=[]

                #Get transcript ids ver
                for l in thisdata:
                    thistx=l.split('\t')[0]
                    if '|' in thistx: thistx=thistx.split('|')[0]
                    txids.append(thistx)

                dftpm=pd.DataFrame({'TranscriptID':txids})
                dfcount=pd.DataFrame({'TranscriptID':txids})


                #current quant file
                for qf in input:
                        # current RAid
                        name=qf.split('/')[1]
                        names.append(name)
                        #Get TPM
                        tpm=pd.read_csv(qf,sep='\t',usecols=[3],skiprows=0)
                        #Get Counts
                        counts=pd.read_csv(qf,sep='\t',usecols=[4],skiprows=0)
                        #Add TPM and counts to respective df with current Run accession ID as column name
                        dftpm[name]=tpm['TPM']
                        dfcount[name]=counts['NumReads']

                ##transcript level TPMs and counts.
                #Read in metadata
                tx_md=pd.read_csv('Transcript_level_metadata.tsv',sep='\t',skiprows=0)
                #Fetch transcript id vers and TPMs for all RA ids.
                df_tpm=dftpm[['TranscriptID']+names].copy()
                df_count=dfcount[['TranscriptID']+names].copy()
                #Merge with metadata
                df_tpm = tx_md.merge(df_tpm, on=['TranscriptID'])
                df_count = tx_md.merge(df_count, on=['TranscriptID'])
                df_tpm.to_csv(output[0],sep='\t',index=False)
                df_count.to_csv(output[1],sep='\t',index=False)


                #Read in gene metadata
                md=pd.read_csv('Gene_level_metadata.tsv',sep='\t',skiprows=0)
                df_gene_tpm=dftpm[['TranscriptID']+names].copy()
                df_gene_count=dfcount[['TranscriptID']+names].copy()

                #Make df with transcript id ver and corresponding Gene_stable_ID ver.
                md2=tx_md[['TranscriptID','Gene_ID_ver']]
                ##Merge TPM and count data with ids
                df_gene_tpm=md2.merge(df_gene_tpm, on=['TranscriptID'])
                df_gene_count=md2.merge(df_gene_count, on=['TranscriptID'])


                #Collapse so that each gene id is listed once. sum up corresponding transcript TPM and counts.
                df_gene_tpm.drop('TranscriptID', axis=1, inplace=True)
                df_gene_count.drop('TranscriptID', axis=1, inplace=True)
                df_gene_tpm = df_gene_tpm.groupby(['Gene_ID_ver'],as_index = False).sum()
                df_gene_count = df_gene_count.groupby(['Gene_ID_ver'],as_index = False).sum()


                ##Merge metadata to counts and TPM
                df_gene_tpm=md.merge(df_gene_tpm, on=['Gene_ID_ver'])
                df_gene_count=md.merge(df_gene_count, on=['Gene_ID_ver'])

                #reorder so gene name is first.
                df_gene_tpm = df_gene_tpm[ ['Gene_name'] + [ col for col in df_gene_tpm.columns if col != 'Gene_name' ] ]
                df_gene_count = df_gene_count[ ['Gene_name'] + [ col for col in df_gene_count.columns if col != 'Gene_name' ] ]



                df_gene_tpm.to_csv(DIR+'/results_TPM_gene.tsv',sep='\t',index=False)
                df_gene_count.to_csv(DIR+'/results_Count_gene.tsv',sep='\t',index=False)
