import pandas as pd
import glob






#Will appened Gene_stable_IDs for each combonation
combos=[]

#Files that contain the EB genes that passed each combonation filter.
files=glob.glob("*_Covid_*")

for i in files:
    df=pd.read_csv(i,sep='\t')
    #Sample # must be >= 5 and be EB
    if (len([col for col in df if col.startswith('batch')])) >=5:
        combos.append(df[((df['Gene_type'] == 'EB_novel'))][["Gene_stable_ID"]])


#Unique list of covid EBs
covidEBs=pd.concat(combos).drop_duplicates()
covidEBs.to_csv("covidEBs.txt",sep='\t',index=False,mode='w')




#Will appened Gene_stable_IDs for each combonation
combos=[]

#Files that contain the EB genes that passed each combonation filter.
files=glob.glob("*_No-Covid_*")

for i in files:
    df=pd.read_csv(i,sep='\t')
    #Sample # must be >= 5 and be EB
    if (len([col for col in df if col.startswith('batch')])) >=5:
        combos.append(df[((df['Gene_type'] == 'EB_novel'))][["Gene_stable_ID"]])


#Unique list of non-covid EBs
nonCovidEBs=pd.concat(combos).drop_duplicates()
nonCovidEBs.to_csv("nonCovidEBs.txt",sep='\t',index=False,mode='w')
