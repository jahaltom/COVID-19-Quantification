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
* SRA_Quant: This code takes in run_accession IDs from the sequencing read archive(place in RAids.txt) and will quantify the expression on Human and SARS-COV-2 transcriptome as well as EB transcripts. Outputs TPM and counts. 

* MultiStudy_Filter.py: This takes output from SRA_Quant(results_TPM_gene.tsv and results_Count_gene.tsv) and separates TPM and counts by study. Additionaly, the TPM is filtered by removing those EB genes where the median < 1, and then removes EB genes where median TPM is less than that of the median of all lncRNA genes. The same genes are removed from counts as well. A tab delimited file called "List" that contans the study_accession IDs and their corrisponding run_accession IDs in that order is also needed. This script also outputs a file that contains all the genes that passed the TPM filter for any study and the corresponding counts as well. 
 
* Bam-Fastq_Quant: Takes in bam/fastq files and will quantify the expression on Human and SARS-COV-2 transcriptome as well as EB transcripts. Outputs TPM and counts.

Bam/Fastq files should be in this structre:
out/Sample1/Sample1.fastq
out/Sample2/Sample2.fastq

* SingleStudy_Filter.py: This takes output from Bam-Fastq_Quant(results_TPM_gene.tsv and results_Count_gene.tsv) and filters it by removing EB genes where the median TPM < 1.  


## Execution 
snakemake -j 50 -s snakefile --cluster "sbatch -t 01:00:00 -c 25"



This code is a modified version of urmi-21's code: https://github.com/urmi-21/covid-quantification.git





