import pandas as pd
import os
import glob

os.chdir(r"C:\Users\15154\Desktop\MesAur\DESeq2")

#Read in DGE-HumanOrtho files
files=glob.glob("*HamsterHumanOrtho*")

#Read in mean of DESeq2 normalized expression file. Round to 3sf
mean_df=pd.read_csv("meanInfo.tsv",sep="\t")
mean_df=mean_df.round(decimals=3)
with pd.ExcelWriter('Syrian_Golden_Hamster_PRJNA837993.xlsx') as writer:
    for i in files:    
        df1=pd.read_csv(i,sep='\t')
        #Rounf every number to 3sf and sort by adjusted p-value.
        df1=df1.round(decimals=3)
        df1=df1.sort_values('padj')
 
        #Gather column names to be extracted from mean_df. 
        cond=i.split("_vs_")[0]+"_Mean"
        condMock=i.split("_vs_")[1].split("_.DGE")[0]+"_Mean"
        #Insert temp column for merge
        df1["ensembl_gene_id_version"]=df1['mauratus_ensembl_gene_id_version'].fillna(df1['SarsCov2_ensembl_gene_id_version'])
        #Subset mean_df   
        df2=mean_df[["ensembl_gene_id_version",cond,condMock]]
        #Merge and filter
        df1=pd.merge(df1,df2,on=["ensembl_gene_id_version"])
    
        df1=df1[['Human_HGNC.symbol', 'Human_Gene.stable.ID', 'mauratus_ensembl_gene_id',
        'mauratus_gene_biotype', 'mauratus_description', cond,condMock,'baseMean',
        'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj','mauratus_ensembl_gene_id_version', 'mauratus_Human.homology.type','Human_Gene.type',
        'Human_Gene.description','SarsCov2_Gene_name','SarsCov2_ensembl_gene_id_version',
       'SarsCov2_ensembl_transcript_id_version',
       'SarsCov2_Gene_description', 'SarsCov2_Gene_type', 'SarsCov2_chr',
       'SarsCov2_seq']]  

    
    
        i=i.replace("-CoV-2_(1000_PFU)", "").replace("_(100,000_PFU)", "").replace("uenza", "").replace("_dpi", "dpi")
        i=i.split("_vs_Mock")[0] 
                
        df1.to_excel(writer, sheet_name=i,index=None)
         
        # Get the xlsxwriter workbook and worksheet objects.
        workbook  = writer.book
        worksheet = writer.sheets[i]
          
        # Add a header format.
        header_format = workbook.add_format({
            'bold': True,
            'valign': 'top',
            'fg_color': '#D7E4BC',
            'border': 1})
          
        # Write the column headers with the defined format.
        for col_num, value in enumerate(df1.columns.values):
            worksheet.write(0, col_num, value, header_format)
            column_length = max(df1[value].astype(str).map(len).max(), len(value))  
            col_idx = df1.columns.get_loc(value)
            worksheet.set_column(col_idx, col_idx, column_length)
            
 
    


