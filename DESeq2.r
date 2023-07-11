library("tximport")
library(DESeq2)


#Make tx2gene file
txMD = read.table("Transcript_level_metadata.tsv",header=TRUE,sep = '\t',quote="")
tx2gene=txMD[,c("TranscriptID","Gene_ID_ver")]




#List quant.sf files
samples = read.table("ids.txt",header=FALSE,sep = '\t')
files <- file.path("out", samples$V1, "salmon_out", "quant.sf")
names(files) <- paste0(samples$V1)

#tximport
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

#Read in metadata
md = read.table("metadata.tsv",header=TRUE,sep = '\t')
rownames(md)=md$SampleID
#Reorder metadata rows to match count data col
md=md[colnames(txi.salmon$counts),]
#Make DeSeq object
dds <- DESeqDataSetFromTximport(txi.salmon,colData= md, design = ~Condition)

##Filter DeSeq object
#filteredGeneList = read.table("results_Count_gene.filtered.tsv",header=TRUE,sep = '\t',quote="",check.names=FALSE)
#Obtain the indices of only desired genes
#genesTokeep <- which(rownames(dds) %in% filteredGeneList$Gene_ID_ver)
#extract your desired genes in the DESeq object
#dds <- dds[genesTokeep, ]

##DESeq2
dds = DESeq(dds)

#Output normalized counts
normCounts=counts(dds, normalized=TRUE)
normCounts = cbind(Gene_ID_ver = rownames(normCounts), normCounts)
write.table(normCounts,"NormalizedCounts.tsv" ,sep = '\t',row.names = FALSE)

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
