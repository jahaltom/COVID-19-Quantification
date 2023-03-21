library("tximport")
library(DESeq2)


#Make tx2gene file
txMD = read.table("Transcript_level_metadata.tsv",header=TRUE,sep = '\t')
tx2gene=data.frame(txMD$TranscriptID,txMD$Gene_ID_ver)

#List quant.sf files
samples = read.table("ids.txt",header=FALSE,sep = '\t')
files <- file.path("out", samples$V1, "salmon_out", "quant.sf")
names(files) <- paste0(samples$V1)

#tximport
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

#Read in metadata
md = read.table("metadata.tsv",header=TRUE,sep = '\t',row.names=2)
#Reorder metadata rows to match count data col
md=md[colnames(txi.salmon$counts),]
#Make DeSeq object
dds <- DESeqDataSetFromTximport(txi.salmon,colData= md, design = ~Condition)

##Filter DeSeq object
filteredGeneList = read.table("results_Count_gene.filtered.tsv",header=TRUE,sep = '\t')
#Obtain the indices of only desired genes
genesTokeep <- which(rownames(dds) %in% filteredGeneList$Gene_ID_ver)
#extract your desired genes in the DESeq object
dds <- dds[genesTokeep, ]

##DESeq2


dds = DESeq(dds)

gene_metadata = read.table("Gene_level_metadata.tsv",header=TRUE,sep = '\t',quote="")


list = read.table("Combo.tsv",sep = '\t',quote="")

for (x in 1:nrow(list)) {
    #########
    #Contrast
    result = results(dds, contrast=c("Condition",list[x,1],list[x,2]))
    ## Remove rows with NA
    result = result[complete.cases(result),]
    #Put GeneID as column
    result = cbind(Gene_ID_ver = rownames(result), result)
    result=merge(gene_metadata,result,by="Gene_ID_ver",all.y=TRUE)
    write.table(result,paste(list[x,1],"vs",list[x,2],".DGE.tsv",sep='_') ,sep = '\t',row.names = FALSE)

}
