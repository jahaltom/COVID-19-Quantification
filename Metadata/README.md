# Making the metadata

Gather transcript_level_metadata.tsv(phylostrata),gene_level_metadata.tsv(phylostrata),EB_Gencode_pep_NODOTS.fasta (Human GencodeV36, EB gene peptide seq), Sars.pep(SARS-COV-2 ASM985889v3 peptide seq), and meta_data.tsv(general metadata for EB,Human annotated,SARs). 

Create fasta files for human
* Run: `conda bash make.fasta.sh` 

Create final fasta
* Run: `python make_metadata.py`


## Metadata
* meta_data.tsv.gz: A metadata file for Human annotated, SARS-COV-2, as well as the EB genes used in this study. 
