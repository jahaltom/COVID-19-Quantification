library(DESeq2)
library(dplyr)

#Read in count information.
countData = read.table("results_Count_gene.filtered.tsv",header=TRUE,sep = '\t',quote="",check.names=FALSE)

#X an Y gene names can be the same. This keeps 1st occurance.
countData <- countData[!duplicated(countData$Gene_stable_ID),]
rownames(countData)= countData$Gene_stable_ID


#Select only columns that contain counts.
countData=select(countData,contains("batch"))

#Round to nearest int
countData=round(countData,0)

##Read in expermental design and gene metadata. rownames is SampleID.
metadata = read.table("metadata.tsv",header=TRUE,row.names=1,sep = '\t')
gene_metadata = read.table("Gene_level_metadata.tsv",header=TRUE,sep = '\t')



#Reorder count columns to match metadata
countData=countData[,rownames(metadata)]


#Should return TRUE
#all(rownames(metadata) == colnames(countData))





##############################EB
##Make DEseq2 object
dds = DESeqDataSetFromMatrix(countData = countData,colData = metadata,design = ~ SequencingBatch + Condition)
dds = DESeq(dds)



list = read.table("Combo.tsv",sep = '\t')

for (x in 1:nrow(list)) {
    #########
    #Contrast
    result = results(dds, contrast=c("Condition",list[x,1],list[x,2]))
    ## Remove rows with NA
    result = result[complete.cases(result),]
    #Put GeneID as column
    result = cbind(Gene_stable_ID = rownames(result), result)
    result=merge(gene_metadata,result,by="Gene_stable_ID",all.y=TRUE)
    write.table(result,paste(list[x,1],"vs",list[x,2],".DGE.tsv",sep='_') ,sep = '\t',row.names = FALSE)

}



