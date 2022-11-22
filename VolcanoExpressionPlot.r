library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)

can=unlist(read.csv('Canonical.txt',sep='\t',header=FALSE))
non_can=unlist(read.csv('Non-Canonical.txt',sep='\t',header=FALSE))
gene_md=read.csv('Gene_level_metadata.tsv',sep='\t')


filelist = list.files(pattern = "DGE.tsv")    
for (file in filelist){
        
    #File name with extension removed
    s=substring(file, 1, nchar(file) - 9)
    #Read in DESeq2 data
    degs=read.csv(file,sep='\t')
    #Merge with metadata
    degs=merge(gene_md,degs,by="Gene_stable_ID",all.y=TRUE)
    #Extract specific immune genes.
    cann=degs[degs$Gene_stable_ID %in% can,]
    non_cann=degs[degs$Gene_stable_ID %in% non_can,]
    immune_genes=rbind(cann,non_cann)
    
    
    #Extract orphan genes
    degs= degs[degs$final_strata.x == "Homo sapiens" & (degs$EB_type == "EB_novel" | degs$EB_type  == "COVID-19 Related EB_novel"),]
    #Concat
    degs=rbind(degs,immune_genes)
    
    

    
    
    # create custom key-value pairs 
    keyvals <- ifelse(
      (degs$EB_type  == "EB_novel"), 'gold',
        ifelse((degs$EB_type  == "COVID-19 Related EB_novel"), 'royalblue',
          'black'))
    keyvals[is.na(keyvals)] <- 'black'
    names(keyvals)[keyvals == 'gold'] <- 'Orphan Gene'
    names(keyvals)[keyvals == 'black'] <- 'Immune Gene'
    names(keyvals)[keyvals == 'royalblue'] <- 'COVID-19 Related Orphan Gene'
    
    
    
    
    
    EnhancedVolcano(degs,
        lab = degs$Gene_name.x,
        x = 'log2FoldChange',
        y = 'padj',
        labCol='black',
        ylim = c(0, 6),
        xlim = c(-10, 10),
        selectLab = unlist(immune_genes$Gene_name.x),
        max.overlaps=100,
        xlab = bquote(~Log[2]~ 'fold change'),
        title = s,
        pCutoff = 0.05,
        FCcutoff = 2,
        pointSize = 2.5,
        colCustom = keyvals,
        colAlpha = 1,
        legendPosition = 'left',
        legendLabSize = 15,
        legendIconSize = 5.0,
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        colConnectors = 'black',
        arrowheads = FALSE,
        labSize=3,
        gridlines.minor = FALSE,
        border = 'partial',
        borderWidth = 1.5,
        borderColour = 'black')
    
    ggsave(paste(s,".Volcano.png"),scale=1,width=20)
    
    
    
    
    
    ggplot(degs, aes(x= baseMean, y = log2FoldChange,colour =EB_type)) + 
    
    geom_point() + xlim(0, 3000) + ylim(-10, 10) + geom_text_repel(data=subset(degs, degs$EB_type  != "EB_novel" & degs$EB_type  != "COVID-19 Related EB_novel" ),
                aes(baseMean,log2FoldChange,label=Gene_name.x), max.overlaps=100) + labs(x="Mean of Normalized Counts") + scale_color_manual(values=c("black", "royalblue", "gold"), 
                           name="Gene Type",
                           labels=c("Immune Gene", "COVID-19 Related Orphan Gene", "Orphan Gene"))
    
    
    ggsave(paste(s,".Expression_L2FC.png"),scale=1,width=20)
    
    
    
    
    
    
    
    
    
}    
    
    
    
