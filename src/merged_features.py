import pandas as pd 
import numpy as np
import sys

lnc_file=sys.argv[1]
prot_file=sys.argv[2]
output_file=sys.argv[3]


prot_df=pd.read_csv(prot_file,sep='\t')


lncrna_df=pd.read_csv(lnc_file,sep='\t')

final=pd.merge(lncrna_df,prot_df,on=['Protein'],how='inner')

final.drop(['Protein'],axis=1,inplace=True)


for cols in final.columns[1:]:
    final[cols]=final[cols].astype(np.float16)

final.to_csv(output_file,index=False,sep='\t')




