# COVID-19-Quantification
This pipeline will take transcriptome data and quantify the expression on the Human (GencodeV36) and SARS-COV-2(ASM985889v3) transcriptomes as well as human EB gene treanscripts using pyrpipe https://github.com/urmi-21/pyrpipe.

## Prerequisites
* Run: `conda env create -f env.yaml`
* Activate the environment: `conda activate pyrpipe_covid`


## Prepareing the reference data. 
* Prepare data: Builds a salmon index using the Human transcriptome(GencodeV36), SARS-COV-2 transcriptome (ASM985889v3), and Human EB gene transcripts. The Human Genome along with viral decoysand spike-insfrom the Genomic Data Commons were also used with Salmon.


## Parameters
The salmon tool parameters are specified in the params directory(needed by pyrpipe). To modify the parameters edit salmon_index.yaml. pyrpipe parameters are specified in the pyrpipe_conf.yaml file. config.yaml contains important parameters for the Snakefiles(SRA_Data and Mason_Data) such as thread number for Salmon, path to fasta, output directory name. 


## Snakefiles and filters
* SRA_Data: This code takes in run_accession IDs from the sequencing read archive(place in RAids.txt) and will quantify the expression on Human and SARS-COV-2 transcriptome as well as EB transcripts. Outputs TPM and counts. 

* SRA_Data_Filter:This takes output from Snakefile(results_TPM_gene.tsv and results_Count_gene.tsv) and separates TPM and counts by study. Additionaly, the TPM is filtered by removing those EB genes where the median < 1, and then removes EB genes where median TPM is less than that of the median of all lncRNA genes. A file called "List" that contans the study_accession IDs and their corrisponding run_accession IDs is also needed.
 
* Mason_Data: Takes bam files from the Mason Nasal Covid study and converts to fastq, then it will quantify the expression on Human and SARS-COV-2 transcriptome as well as EB transcripts. Outputs TPM and counts.

* Mason_Data_Filter: This takes output from Mason Data(results_TPM_gene.tsv) and filters teh TPM by removing those protein coding and EB genes where the median < 1, and then removes EB genes where median TPM is less than that of the median of all lncRNA genes. 


## Execution 
snakemake -j 50 -s snakefile --cluster "sbatch -t 01:00:00 -c 60"



This code is a modified version of urmi-21's code: https://github.com/urmi-21/covid-quantification.git





