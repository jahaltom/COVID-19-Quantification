# COVID-19-Quantification
This pipeline will take Human RNA-Seq data and concurrently quantify the expression on the Human (GencodeV36) and SARS-COV-2(ASM985889v3) transcriptomes as well as human evidence based (EB) gene treanscripts using pyrpipe https://github.com/urmi-21/pyrpipe.

## Prerequisites
* Run: `conda env create -f env.yaml`
* Activate the environment: `conda activate pyrpipe_covid`


## Prepareing the reference data. 
* Prepare data: Builds a salmon index using the Human transcriptome(GencodeV36), SARS-COV-2 transcriptome (ASM985889v3), and Human EB gene transcripts. To increase mapping accuracy, the human genome along with viral decoys and spike-ins from the Genomic Data Commons (GRCh38.d1.vd1) will be used as decoys for Salmon.


## Parameters
pyrpipe parameters are specified in the pyrpipe_conf.yaml file. config.yaml contains important parameters for the Snakefiles(Bam-Fastq_Quant.py and SRA_Quant.py ) such as thread number for Salmon, path to fasta/bam, output directory name. 


## Snakefiles, QC, and filters
* SRA_Quant.py: This code takes in a single column list of run_accession IDs from the sequencing read archive(place in ids.txt) and will quantify the expression on Human and SARS-COV-2 transcriptome as well as EB transcripts. Outputs TPM and counts at gene and transcript level. Gene level summed up across transcripts. Run FastQC in the default mode on the FastQ file once trimming is complete

* Bam-Fastq_Quant.py: Takes in bam/fastq files and will quantify the expression on Human and SARS-COV-2 transcriptome as well as EB transcripts. Outputs TPM and counts at gene and transcript level. Gene level summed up across transcripts. Run FastQC in the default mode on the FastQ file once trimming is complete

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

* MultiQC
```
mkdir multiqc

#Move fastqc zip files into single directory
cat ids.txt | while read i;do 
	mv out/$i/*zip* multiqc/; 
done

#Run MultiQC
multiqc multiqc/*

#Same for salmon output 

mkdir multiqc/salmon

cat ids.txt | while read i;do 
	cp -r out/$i/*salmon_out* multiqc/salmon/salmon_out$i; 
done

multiqc -n multiqc_report_rnaseq multiqc/salmon/
```

* Filter.py: This takes output from SRA_Quant.py or Bam-Fastq_Quant.py (results_TPM_gene.tsv and results_Count_gene.tsv) and separates TPM and counts into groups by an experimental condition of your choosing. Per group, the TPM is then filtered by removing those genes where the median < 1 (except SARS-COV-2 transcripts). EB genes are removed where the median TPM is less than that of the median of the medians of all known protein coding genes. EB genes from Urmis 54K list are excluded/or not from filtering in (See hastag at #Remove EB genes where median TPM is less than that of the med_med_pc). The same genes are removed from counts as well. A tab delimited file called "metadata.tsv" that contains the Condition and SampleID is also needed. This script also outputs a file that contains the TPM for all the genes that passed the TPM filter for any condition and the corresponding counts as well (results_TPM_gene.filtered.tsv results_Count_gene.filtered.tsv). For each condition, outputs the median of the medians of the known protein coding, lincRNA, and EB genes (medianInfo.tsv) and creates a TPM tsv of genes that made it through each condition in the filter *_EB.tsv*. 
 

* SeparateByStudy.py: This takes output from SRA_Quant.py (results_TPM_gene.tsv and results_Count_gene.tsv). Subsets expression matrixs into their respective studies. Needs a Study.tsv file that contains study_accession and run_accession IDs.
 

* NewEBs.py: Takes in _EB.tsv files from Filter.py for for each Condition(metadata.tsv) that has a sample size >=5 and outputs a list of all covid and non-covid EB genes (Gene_stable_ID) (nonCovidEBs.txt and covidEBs.txt). ##Needs work


## Execution
```
snakemake -j 50 -s snakefile --cluster "sbatch -t 01:00:00 -c 25"
```

## DESeq2
* DESeq2.r:
Takes in Transcript/Gene _level_metadata.tsv, ids.txt, quant.sf files (from snakefiles), metadata.tsv, and Combo.tsv.
Combo.tsv: is a file that is 2 column and tab spaced(Condition1	Condition2). For each line the constrast is made from left to right (Condition1	vs Condition2).
Can supply expression matrix (results_Count_gene.filtered.tsv) to filter DeSeq object.


## PCA
needs work


