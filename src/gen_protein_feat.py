
import pandas as pd 
import numpy as np
import itertools
from tqdm import tqdm
import re
import sys

prot_name=sys.argv[1]
output_file=sys.argv[2]




# Amino acides - They are only 20 
amino_acides = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']



def prot_feat(prot_df,props,segre_list):
    li=segre_list
    li_2 = []
    for i in itertools.product(segre_list, repeat=2):
        li_2.append(''.join(map(str, i)))

    li_3 = []
    for i in itertools.product(segre_list, repeat=3):
        li_3.append(''.join(map(str, i)))

    
    
    for i in tqdm(range(0,prot_df.shape[0])):
        for j in li:
            prot_df.loc[i,j+'_'+props+'_prot_fraction']=prot_df.loc[i,'Protein_Sequence'].count(j)/len(prot_df.loc[i,'Protein_Sequence'])



        for j in li_2:
            prot_df.loc[i,j+'_'+props+'_prot_fraction']=len(re.findall('(?='+j+')',prot_df.loc[i,'Protein_Sequence']))/(len(prot_df.loc[i,'Protein_Sequence'])-1)
            


        for j in li_3:
            prot_df.loc[i,j+'_'+props+'_prot_fraction']=len(re.findall('(?='+j+')',prot_df.loc[i,'Protein_Sequence']))/(len(prot_df.loc[i,'Protein_Sequence'])-2)

            
    return prot_df
    



df=pd.read_csv('../dataset/Final_Protein.txt',sep='\t')
df.drop(['Prot_ID'],axis=1,inplace=True)
df=df[df.Protein==prot_name]
df['Protein_Sequence']=df['Protein_Sequence'].str.upper()
df.reset_index(inplace=True,drop=True)


hydrophobic = ['A', 'C', 'I', 'L', 'M', 'F', 'W', 'V']
neutral = ['G', 'H', 'P', 'S', 'T', 'Y']
hydrophilic = ['R', 'N', 'D', 'Q', 'E', 'K']

segre_list=['1','2','3']

hydropathy_df=df.copy()


for amn in hydrophobic:
    hydropathy_df['Protein_Sequence']=hydropathy_df.Protein_Sequence.str.replace(amn,'1')
for amn in neutral:
    hydropathy_df['Protein_Sequence']=hydropathy_df.Protein_Sequence.str.replace(amn,'2')
for amn in hydrophilic:
    hydropathy_df['Protein_Sequence']=hydropathy_df.Protein_Sequence.str.replace(amn,'3')

hydropathy_df=prot_feat(hydropathy_df,'hydropathy',segre_list)
hydropathy_df.drop(['Protein_Sequence'],inplace=True,axis=1)




aliphatic=['G', 'A', 'V', 'L', 'I', 'P']
aromatic= ['F', 'Y', 'W']
polar_neutral= ['S', 'T', 'N', 'Q']
sulpher_aa = ['C', 'M']
pos_neg_charged= ['D', 'E', 'H', 'K', 'R']

segre_list=['1','2','3','4','5']

chemical_df=df.copy()


for amn in aliphatic:
    chemical_df['Protein_Sequence']=chemical_df.Protein_Sequence.str.replace(amn,'1')
for amn in aromatic:
    chemical_df['Protein_Sequence']=chemical_df.Protein_Sequence.str.replace(amn,'2')
for amn in polar_neutral:
    chemical_df['Protein_Sequence']=chemical_df.Protein_Sequence.str.replace(amn,'3')
for amn in sulpher_aa:
    chemical_df['Protein_Sequence']=chemical_df.Protein_Sequence.str.replace(amn,'4')
for amn in pos_neg_charged:
    chemical_df['Protein_Sequence']=chemical_df.Protein_Sequence.str.replace(amn,'5')
    
    
chemical_df=prot_feat(chemical_df,'chemical',segre_list)
chemical_df.drop(['Protein_Sequence'],inplace=True,axis=1)



very_small=['A', 'G', 'S']
small=['N', 'D', 'C', 'P', 'T']
medium=['Q', 'E', 'H', 'V']
large=['R', 'I', 'L', 'K', 'M']
very_large=['F', 'W', 'Y']

segre_list=['1','2','3','4','5']

volume_df=df.copy()

for amn in very_small:
    volume_df['Protein_Sequence']=volume_df.Protein_Sequence.str.replace(amn,'1')
for amn in small:
    volume_df['Protein_Sequence']=volume_df.Protein_Sequence.str.replace(amn,'2')
for amn in medium:
    volume_df['Protein_Sequence']=volume_df.Protein_Sequence.str.replace(amn,'3')
for amn in large:
    volume_df['Protein_Sequence']=volume_df.Protein_Sequence.str.replace(amn,'4')
for amn in very_large:
    volume_df['Protein_Sequence']=volume_df.Protein_Sequence.str.replace(amn,'5')


    
volume_df=prot_feat(volume_df,'volume',segre_list)
volume_df.drop(['Protein_Sequence'],inplace=True,axis=1)




gen_df=df.copy()
for aa in amino_acides:
    gen_df[aa+'_fraction']=gen_df['Protein_Sequence'].str.count(aa)/gen_df['Protein_Sequence'].str.len()



gen_df.drop(['Protein_Sequence'],inplace=True,axis=1)


Final=pd.merge(gen_df,hydropathy_df,how='inner',on='Protein')
Final=pd.merge(Final,chemical_df,how='inner',on='Protein')
Final=pd.merge(Final,volume_df,how='inner',on='Protein')


Final.to_csv(output_file,index=False,sep='\t')





