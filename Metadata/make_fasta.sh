#!/bin/bash

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EB_Gencode_pep_NODOTS.fasta | sed 'N;s/\n/ /' | perl -p -i -e 's/ /\t/g' | grep "|ENSG0" > human
cat human | sed 's/|/ /g' | awk 'BEGIN {OFS = "\t"}{print $2,$NF}' > Human
rm human
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EB_Gencode_pep_NODOTS.fasta | sed 'N;s/\n/ /' | perl -p -i -e 's/ /\t/g' | grep -v "|ENSG0" | sed 's/>mikado/EB/g' | awk 'BEGIN {OFS = "\t"}{print $1,$3}' > EB

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Sars.pep | sed 'N;s/\n/ /' |  perl -p -i -e 's/ /\t/g' | sed 's/>//g' > sars