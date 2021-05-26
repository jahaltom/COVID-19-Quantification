# COVID-19-Quantification
This pipeline will take transcriptome data and quantify the expression on Human transcriptome as well as SARs and Human EB genes. 

##Prerequisites
* Run: `conda env create -f env.yaml`
* Activate the environment: `conda activate pyrpipe_covid`


Prepare data at the moment will create a fasta file of the human and SARs transcriptome, and includes Human genome as decoy for salmon. Will use Human genome to make decoys.txt file for salmon.--Can add yur own EB genes. Can also add viral decoys and spike-ins as decoys. 

Snakefile: This code takes in run_accession IDs from the sequencing read archive and will quantify the expression on Human transcriptome as well as SARs and Human EB genes. Outputs TPM and counts. 
Snakefile Filter:This takes output from Snakefile(results_TPM_gene.tsv and results_Count_gene.tsv) and separates TPM and counts by study. Additionaly, the TPM is filtered by removing those protein coding and EB genes where the median < 1, and then removes EB genes where median TPM is less than that of the median of all lncRNA genes. 

Mason Data: Takes bam files from Mason Covid study and converts to fastq, then it will quantify the expression on Human transcriptome as well as SARs and Human EB genes. Outputs TPM and counts.
Mason Data Filter: This takes output from Mason Data(results_TPM_gene.tsv) and filters teh TPM by removing those protein coding and EB genes where the median < 1, and then removes EB genes where median TPM is less than that of the median of all lncRNA genes. 

meta_data.tsv.gz: A metadata file for Human and SARs transcriptome, as well as the EB genes used in this study. 



This code is a modified version of urmi-21's code: https://github.com/urmi-21/covid-quantification.git





