import pandas as pd 
import numpy as np
import sys
import os
import shutil
import matplotlib.pyplot as plt
import logomaker as lm

input_file=sys.argv[1]
img_file=sys.argv[2]
window_len=int(sys.argv[3])


df=pd.read_csv(input_file,sep='\t')
df1=df.copy()
df=df[df.Model_Score>=0.5]


if(df.shape[0]>1):
    df=df[df.lncRNA_Sequence.str.len()==window_len]
    

if (df.shape[0]>0):
    counts_mat = lm.alignment_to_matrix(df.lncRNA_Sequence.values)
    logo_sav=lm.Logo(counts_mat)
    prot_name=df.Protein.unique()[0]
    inter_num=str(np.round((df.shape[0]*100/df1.shape[0]),2))
    # title='Motif plot corresponding to '+prot_name+' protein containing '+inter_num+'% positive interaction'
    title='Sequence Logo representing Consensus Sequence Motif within the '+prot_name+' Protein specific lncRNA Interacting Segments'
    plt.title(title)
    plt.savefig(img_file)

else:
    src_file='../img/no_interaction.png'
    shutil.copy(src_file,'../output/')

    if os.path.exists(img_file):
        os.remove(img_file)
    os.rename('../output/no_interaction.png',img_file)




