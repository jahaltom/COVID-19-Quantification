# COVID-19-Quantification
This pipeline will take Human RNA-Seq data and concurrently quantify the expression on the Human (GencodeV36) and SARS-COV-2(ASM985889v3) transcriptomes as well as human evidence based (EB) gene treanscripts using pyrpipe https://github.com/urmi-21/pyrpipe.

## Prerequisites
* Run: `conda env create -f env.yaml`
* Activate the environment: `conda activate pyrpipe_covid`


## Prepareing the reference data. 
* Prepare data: Builds a salmon index using the Human transcriptome(GencodeV36), SARS-COV-2 transcriptome (ASM985889v3), and Human EB gene transcripts. To increase mapping accuracy, the human genome along with viral decoys and spike-ins from the Genomic Data Commons (GRCh38.d1.vd1) will be used as decoys for Salmon.


## Parameters
The salmon tool parameters are specified in the params directory(needed by pyrpipe). To modify the parameters edit salmon_index.yaml. pyrpipe parameters are specified in the pyrpipe_conf.yaml file. config.yaml contains important parameters for the Snakefiles(SRA_Data and Mason_Data) such as thread number for Salmon, path to fasta, output directory name. 


## Snakefiles and filters
* SRA_Quant.py: This code takes in a single column list of run_accession IDs from the sequencing read archive(place in RAids.txt) and will quantify the expression on Human and SARS-COV-2 transcriptome as well as EB transcripts. Outputs TPM and counts at gene and transcript level. Gene level summed up across transcripts. 

* Filter.py: This takes output from SRA_Quant.py or Bam-Fastq_Quant.py (results_TPM_gene.tsv and results_Count_gene.tsv) and separates TPM and counts into groups by an experimental condition of your choosing. Per group, the TPM is then filtered by removing those genes where the median < 1 (except SARS-COV-2 transcripts).  EB genes are removed where the median TPM is less than that of the median of the medians of all known protein coding genes. EB genes from Urmis 54K list are excluded from filtering. The same genes are removed from counts as well. A tab delimited file called "List.tsv" that contans the Condition and SampleID in that order is also needed. This script also outputs a file that contains the TPM for all the genes that passed the TPM filter for any condition and the corresponding counts as well. Outputs median of known protein coding and lincRNA, and EB genes. Creates a quantification TSV of genes that made it through filter for each condition. 

* SeparateByStudy.py: Subsets an expression matrix into their respective studies. Needs a Study.tsv file that contains study_accession and run_accession IDs in that order.
 
* Bam-Fastq_Quant.py: Takes in bam/fastq files and will quantify the expression on Human and SARS-COV-2 transcriptome as well as EB transcripts. Outputs TPM and counts at gene and transcript level. Gene level summed up across transcripts. 

  * Specify FileType and Layout in config.yaml.
  * ids.txt must contain directory names for individual sample fastq(s). Below would be (Sample1,Sample2) <-one per line in ids.txt. 
  * Bam/Fastq files should be in this structre:
```
#Single-end
out/Sample1/Sample1.fastq
out/Sample2/Sample2.fastq

#Paired-end
out/Sample1/Sample1.r1.fastq
out/Sample1/Sample1.r2.fastq
```



## Execution
```
snakemake -j 50 -s snakefile --cluster "sbatch -t 01:00:00 -c 25"
```






