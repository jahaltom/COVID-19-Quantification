# COVID-19-Quantification
This pipeline will take transcriptome data and quantify the expression on Human transcriptome as well as SARs using pyrpipe https://github.com/urmi-21/pyrpipe.

## Prerequisites
* Run: `conda env create -f env.yaml`
* Activate the environment: `conda activate pyrpipe_covid`


## Prepareing the reference data. 
* Prepare data: This will create a fasta file of the Human and SARs transcriptome, and includes Human genome as decoy for salmon. Will use Human genome to make decoys.txt file for salmon.(For this analysis it also included our Human EB gene trancritpome as well as viral decoys and spike-ins). 


## Parameters
The salmon tool parameters are specified in the params directory(needed by pyrpipe). To modify the parameters edit salmon_index.yaml. pyrpipe parameters are specified in the pyrpipe_conf.yaml file. config.yaml contains important parameters for the Snakefiles(SRA_Data and Mason_Data). 


## Snakefiles and filters
* SRA_Data: This code takes in run_accession IDs from the sequencing read archive(place in RAids.txt) and will quantify the expression on Human and SARs transcriptome. Outputs TPM and counts. 

* SRA_Data_Filter:This takes output from Snakefile(results_TPM_gene.tsv and results_Count_gene.tsv) and separates TPM and counts by study. Additionaly, the TPM is filtered by removing those protein coding and EB genes where the median < 1, and then removes EB genes where median TPM is less than that of the median of all lncRNA genes. A file called "List" that contans the study_accession IDs and their corrisponding run_accession IDs is alos needed.
 
* Mason_Data: Takes bam files from Mason Covid study and converts to fastq, then it will quantify the expression on Human and SARs transcriptome. Outputs TPM and counts.

* Mason_Data_Filter: This takes output from Mason Data(results_TPM_gene.tsv) and filters teh TPM by removing those protein coding and EB genes where the median < 1, and then removes EB genes where median TPM is less than that of the median of all lncRNA genes. 


## Metadata
* meta_data.tsv.gz: A metadata file for Human and SARs genes, as well as the EB genes used in this study. 


## Execution 
snakemake -j 50 -s snakefile --cluster "sbatch -t 01:00:00 -c 60"



This code is a modified version of urmi-21's code: https://github.com/urmi-21/covid-quantification.git





