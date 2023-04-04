library(dplyr)
library(tibble)
library("ggplot2")
library("ggrepel")
library("tximport")
library(DESeq2)

#Make tx2gene file
txMD = read.table("Transcript_level_metadata.tsv",header=TRUE,sep = '\t',quote="")
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





##Filter DeSeq object for common (between Robs and Mine)
commGenes = read.table("IDsToUse",header=FALSE,sep = '\t',quote="",check.names=FALSE)
colnames(commGenes)="Gene_ID_ver"
#Obtain the indices of only desired genes
genesTokeep <- which(rownames(dds) %in% commGenes$Gene_ID_ver)
#extract your desired genes in the DESeq object
ddsComm <- dds[genesTokeep, ]
#Transform
vsd <- vst(ddsComm, blind=FALSE)
rld <- rlog(ddsComm, blind=FALSE)
write.table(assay(vsd),"vsdComm.tsv",sep = '\t',row.names = TRUE)
write.table(assay(rld),"rldComm.tsv",sep = '\t',row.names = TRUE)

library(DESeq2)
library("ggplot2")


vsd = read.table("vsdComm.tsv",header=TRUE,sep = '\t',check.names=FALSE)
md = read.table("metadata.tsv",header=TRUE,sep = '\t',row.names=2)
#Reorder metadata rows to match count data col
md=md[colnames(vsd),]
vsd=round(vsd)
vsd <- DESeqDataSetFromMatrix(vsd,colData= md, design = ~Condition)


plotPCA(vsd, intgroup=c("OE", "Cell"))
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#####################################

#Read in TPM information.
df = read.table("results_TPM_gene.Comm.tsv",header=TRUE,sep = '\t',quote="",row.names=1,check.names=FALSE)
df=log2(df+ 1)

data.pca = prcomp(t(df),center=TRUE)

##Make Run column
data.pca = tibble::rownames_to_column(as.data.frame(data.pca$x), "SampleID")

##Merge with metadata
data.pca.metadata=merge(metadata,data.pca,all.by="SampleID",all.x = TRUE,all.y=TRUE)

ggplot(data.pca.metadata, aes(x = PC1, y = PC2, color = Cell ,shape=OE,label = SampleID)) +
  geom_point(size =2) +
  ggtitle("PC1 vs PC2 Common")  + geom_text_repel(aes(label = SampleID), size = 3.5)
ggsave("PCA.TissuevsInfection.png")

