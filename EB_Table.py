import pandas as pd
from pandas import DataFrame
import glob

#Will appened info for each combonation
combos=[]
#Read in median info.
median=pd.read_csv('medianInfo.tsv',sep='\t',header=None)  
#Files that contain the genes that passed each combonation filter. 
files=glob.glob("*EB*")

for i in files:
    df=pd.read_csv(i,sep='\t')     
    #Combonation, Sample #, Jeff EB Total, Urmi EB Total, Jeff EB Orphan, Urmi EB Orphan. 
    combos.append([i.replace('_EB.tsv' , ''),(len([col for col in df if col.startswith('SRR')]) or len([col for col in df if col.startswith('ERR')])),df['is54K_EB'].value_counts()[0],df['is54K_EB'].value_counts()[1] ,len(df[(df['ps'] == "27.0") & (df['is54K_EB'] == False)]),len(df[(df['ps'] == "27.0") & (df['is54K_EB'] == True)])])
    

combos=DataFrame(combos)
combos.columns =['Combonation', 'Sample #', 'Jeff EB Total', 'Urmi EB Total', 'Jeff EB Orphan','Urmi EB Orphan']
#Combonation and Protein Coding, lincRNA, and Evidence based median of medians
median=median[[0,2,4,6]]
median.columns =['Combonation', 'Protein Coding median of medians', 'lincRNA median of medians', ' Evidence based median of medians']
#Merge
combos=pd.merge(combos,median,on=['Combonation'])


combos.to_csv("combonations.tsv",sep='\t',index=False,mode='w')
