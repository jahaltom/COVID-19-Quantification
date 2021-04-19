##Extract EB genes
cat human_tr_gen_decoy.fasta | grep ">" | grep "gene=EB." |  awk 'BEGIN{OFS="\t"}{print $1,$2,"NA",$NF,"NA","NA","EB"}'  | sed 's/>//g' | sed 's/gene=//g' > EB
rm human_tr_gen_decoy.fasta

Get sars metsdata:
cat Ens_gene_metadata.txt  | grep "ENSSAST" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$5,$1,$2,$3,$4,$7 }' > Sars
rm Ens_gene_metadata.txt

cat gencode_tx_metadata.tsv | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$17,$19,$20,$4,$7,$18}' > head


















#download covid seq
wget -q ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/cdna/Sars_cov_2.ASM985889v3.cdna.all.fa.gz
