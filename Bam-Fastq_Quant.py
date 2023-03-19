import yaml
import sys
import os
import glob
from pyrpipe import sra,quant,qc
from pyrpipe.runnable import Runnable
import pandas as pd


#Read in config file
configfile: "config.yaml"
DIR = config['DIR']
THREADS=config['THREADS']
FileType=config['FileType']
Layout=config['Layout']
Tr=config['Transcriptome']

#creates a trim_galore object.
trim_galore=qc.Trimgalore(threads=4)

#Read in run_accession ids#
with open ("ids.txt") as f:
        ids=f.read().splitlines()

#create salmon object. If index is not present it will be created
salmon=quant.Salmon(index=Tr.split("/")[0]+"/salmon_index",transcriptome=Tr,threads=THREADS)

rule all:
        input:
                expand("{wd}/results_TPM_tx.tsv",wd=DIR)
rule quant:
        output:
                "{wd}/{sample}/salmon_out/quant.sf"
        run:
            
                #Pathway for fastq.gz/fq.gz/bam.gz files and salmon output 
                path=DIR + "/" + wildcards.sample+ "/"
                            
                if (FileType == 'fastq.gz' or FileType == 'fq.gz'):
                    
                    #fastq file paths
                    my_files_path = glob.glob(path+'*.'+FileType)

                    if len(my_files_path) == 1:
                        #Run Salmon on sra object(fastq file) and delete fastq when finished.
                        sra.SRA(fastq=my_files_path[0],directory=path).trim(trim_galore).quant(salmon).delete_fastq()
        
                    elif len(my_files_path) == 2:
                        #Run Salmon on sra object(fastq files) and delete fastqs when finished.
                        sra.SRA(fastq=my_files_path[0],fastq2=my_files_path[1],directory=path).trim(trim_galore).quant(salmon).delete_fastq()                                                               
                        
                elif (FileType == 'bam.gz' and Layout == "Paired"):
                    #bam file paths
                    my_files_path = glob.glob(path+'*.'+FileType)
                    #Fastq files. 
                    fq1=path + wildcards.sample + '1.fastq'
                    #Fastq files. 
                    fq2=path + wildcards.sample + '2.fastq'
                                
                    #Run picard to convert bam to fastq
                    picard=Runnable(command='picard')
                    param={'-Xmx20g':'','SamToFastq':'','I=': my_files_path[0] ,'F=': fq1,'F2=': fq2}
                    picard.run(**param)
                    
                    #Run Salmon on sra object(fastq files) and delete fastqs when finished.
                    sra.SRA(fastq=fq1,fastq2=fq2,directory=path).trim(trim_galore).quant(salmon).delete_fastq()
                elif (FileType == 'bam' and Layout == "Single"):
                    #bam file paths
                    my_files_path = glob.glob(path+'*.'+FileType)
                    #Fastq files. 
                    fq=path + wildcards.sample + '.fastq'
 
                                
                    #Run picard to convert bam to fastq
                    picard=Runnable(command='picard')
                    param={'-Xmx20g':'','SamToFastq':'','I=': my_files_path[0] ,'F=': fq}
                    picard.run(**param)
                    
                    #Run Salmon on sra object(fastq files) and delete fastq when finished.
                    sra.SRA(fastq=fq,directory=path).trim(trim_galore).quant(salmon).delete_fastq()

                    
                    
       
rule merge:
        input:
                ["{wd}/{sample}/salmon_out/quant.sf".format(wd=DIR,sample=s) for s in ids]
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

                #Get transcript ids
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
                #Fetch transcript ids and TPMs for all RA ids.
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

                #Make df with transcript id and corresponding Gene_stable_ID
                md2=tx_md[['TranscriptID','Gene_stable_ID']]
                ##Merge TPM and count data with ids
                df_gene_tpm=md2.merge(df_gene_tpm, on=['TranscriptID'], how='right')
                df_gene_count=md2.merge(df_gene_count, on=['TranscriptID'], how='right')


                #Collapse so that each gene id is listed once. sum up corresponding transcript TPM and counts.
                df_gene_tpm.drop('TranscriptID', axis=1, inplace=True)
                df_gene_count.drop('TranscriptID', axis=1, inplace=True)
                df_gene_tpm = df_gene_tpm.groupby(['Gene_stable_ID'],as_index = False).sum()
                df_gene_count = df_gene_count.groupby(['Gene_stable_ID'],as_index = False).sum()


                ##Merge metadata to counts and TPM
                df_gene_tpm=md.merge(df_gene_tpm, on=['Gene_stable_ID'], how='right')
                df_gene_count=md.merge(df_gene_count, on=['Gene_stable_ID'], how='right')

                #reorder so gene name is first.
                df_gene_tpm = df_gene_tpm[ ['Gene_name'] + [ col for col in df_gene_tpm.columns if col != 'Gene_name' ] ]
                df_gene_count = df_gene_count[ ['Gene_name'] + [ col for col in df_gene_count.columns if col != 'Gene_name' ] ]



                df_gene_tpm.to_csv(DIR+'/results_TPM_gene.tsv',sep='\t',index=False)
                df_gene_count.to_csv(DIR+'/results_Count_gene.tsv',sep='\t',index=False)
