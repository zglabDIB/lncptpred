
import pandas as pd 
import numpy as np
import pickle
import gc
import warnings
import sys

warnings.filterwarnings("ignore")


init_file=sys.argv[1]
feat_file=sys.argv[2]
output_file=sys.argv[3]


scaler=pickle.load(open('../models/StandardScaler_lncRNA_Protein.sav', 'rb'))
model=pickle.load(open('../models/LightGBM_lncRNA_Protein.sav', 'rb'))

init_df=pd.read_csv(init_file,sep='\t')

inp_df=pd.read_csv(feat_file,sep='\t',dtype = 'float16', converters = {'temp_ind': str})

inp_np=inp_df.values
X_inp=inp_np[:,1:]

X_inp=scaler.transform(X_inp)
X_inp=X_inp.astype(np.float16)


model_score=model.predict_proba(X_inp)


inp_df['Model_Score']=np.round(model_score[:,1],3)

inp_df=inp_df[['temp_ind','Model_Score']]


final_df=pd.merge(init_df,inp_df,on=['temp_ind'],how='inner')
final_df.drop(['temp_ind'],axis=1,inplace=True)


final_df.to_csv(output_file,index=False,sep='\t')

