import pandas as pd
import numpy as np
import sys
from tqdm import tqdm


lncrna_file=sys.argv[1]

handle=open(lncrna_file,'r')
for a in handle.readlines():
    lncrna=a

handle.close()





window_len=int(sys.argv[2])
shift_size=int(sys.argv[3])
strand=sys.argv[4]
prot=sys.argv[5]
file_name=sys.argv[6]
total_len=len(lncrna)

df=pd.DataFrame(columns=['Window_Start','Window_End','lncRNA_Sequence','Strand','Protein'])

iteration=int(np.ceil((total_len-window_len)/shift_size))+1

for i in tqdm(range(iteration)):
    start=i*shift_size
    end=i*shift_size+window_len
    df.loc[i,'Window_Start']=start+1
    df.loc[i,'Window_End']=end+1
    df.loc[i,'lncRNA_Sequence']=lncrna[start:end]


df['Strand']=strand
df['Protein']=prot
df['lncRNA_Sequence']=df['lncRNA_Sequence'].str.upper()
df['lncRNA_Sequence']=df['lncRNA_Sequence'].str.replace('U','T')

df['temp_ind']=df['Window_Start'].astype(np.str)+'_'+df['Window_End'].astype(np.str)

df.to_csv(file_name,sep='\t',index=False)




