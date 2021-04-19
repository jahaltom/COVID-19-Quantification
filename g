
Get sars metsdata:
cat Ens_gene_metadata.txt  | grep "ENSSAST" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1".1",$5,$1,$2,$3,$4,$7 }' > Sars
rm Ens_gene_metadata.txt

EB genes and gencode transcripts
cat EB_Gencode_gene_metadata_ordered_withPS.tsv | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$16,$17,$20,$3,$6,$18}' > head
rm EB_Gencode_gene_metadata_ordered_withPS.tsv


cat head Sars  > meta_data.tsv
rm  Sars head

















#download covid seq
wget -q ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/cdna/Sars_cov_2.ASM985889v3.cdna.all.fa.gz
