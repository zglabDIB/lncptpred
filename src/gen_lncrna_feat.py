import pandas as pd 
import numpy as np
import itertools
from tqdm import tqdm
import re
import sys


input_file=sys.argv[1]
output_file=sys.argv[2]



df=pd.read_csv(input_file,sep='\t')
df.loc[df['Strand']=='-','Strand']=-1
df.loc[df['Strand']=='+','Strand']=1
df['Strand']=df['Strand'].astype(np.int8)

li=['A','C','G','T']



li_2 = []
for i in itertools.product(['A','C','G','T'], repeat=2):
    li_2.append(''.join(map(str, i)))


li_3 = []
for i in itertools.product(['A','C','G','T'], repeat=3):
    li_3.append(''.join(map(str, i)))


li_4 = []
for i in itertools.product(['A','C','G','T'], repeat=4):
    li_4.append(''.join(map(str, i)))


for j in tqdm(li):
    df[j+'_lncrna_fraction']=df['lncRNA_Sequence'].str.count(j)/df['lncRNA_Sequence'].str.len()



for j in tqdm(li_2):
    df[j+'_lncrna_fraction']=df['lncRNA_Sequence'].apply(lambda x:len(re.findall('(?='+j+')',x)))/(df['lncRNA_Sequence'].str.len()-1)



for j in tqdm(li_3):
    df[j+'_lncrna_fraction']=df['lncRNA_Sequence'].apply(lambda x:len(re.findall('(?='+j+')',x)))/(df['lncRNA_Sequence'].str.len()-2)


for j in tqdm(li_4):
    df[j+'_lncrna_fraction']=df['lncRNA_Sequence'].apply(lambda x:len(re.findall('(?='+j+')',x)))/(df['lncRNA_Sequence'].str.len()-3)



df.drop(['Window_Start','Window_End','lncRNA_Sequence'],inplace=True,axis=1)


cols_arrange=list(df.columns[1:])+[df.columns[0]]
df=df[cols_arrange]


df.to_csv(output_file,index=False,sep='\t')



