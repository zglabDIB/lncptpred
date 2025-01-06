
import pandas as pd 
import numpy as np
import pickle
from sklearn.preprocessing import StandardScaler
import gc
from sklearn.metrics import accuracy_score
import warnings
from sklearn.metrics import confusion_matrix,classification_report,f1_score,roc_auc_score,auc,roc_curve,matthews_corrcoef
import sys

warnings.filterwarnings("ignore")


init_file=sys.argv[1]
feat_file=sys.argv[2]
output_file=sys.argv[3]


model1=pickle.load(open('../models/LR_Meta_lncRNA_Prot_Balanced.sav', 'rb'))
model2=pickle.load(open('../models/LR_Meta_lncRNA_Prot_Pos_Bias.sav', 'rb'))
model4=pickle.load(open('../models/LR_Meta_lncRNA_Prot_GLOBAL.sav', 'rb'))

scaler1=pickle.load(open('../models/StandardScaler_lncRNA_Protein.sav', 'rb'))

init_df=pd.read_csv(init_file,sep='\t')

inp_df=pd.read_csv(feat_file,sep='\t',dtype = 'float16', converters = {'temp_ind': str})

inp_np=inp_df.values
X_inp=inp_np[:,1:]

X_inp=scaler1.transform(X_inp)
X_inp=X_inp.astype(np.float16)


balanced_pred=model1.predict_proba(X_inp)
pos_biased_pred=model2.predict_proba(X_inp)
averaging_pred=np.mean([balanced_pred,pos_biased_pred],axis=0)
global_pred=model4.predict_proba(X_inp)


inp_df['Balanced_Classifier_Predict_Proba']=np.round(balanced_pred[:,1],3)
inp_df['Balanced_Classifier_Predict']=np.argmax(balanced_pred,axis=1)

inp_df['Positively_Tuned_Classifier_Predict_Proba']=np.round(pos_biased_pred[:,1],3)
inp_df['Positively_Tuned_Classifier_Predict']=np.argmax(pos_biased_pred,axis=1)

inp_df['Averaging_Classifier_Predict_Proba']=np.round(averaging_pred[:,1],3)
inp_df['Averaging_Classifier_Predict']=np.argmax(averaging_pred,axis=1)

inp_df['Global_Classifier_Predict_Proba']=np.round(global_pred[:,1],3)
inp_df['Global_Classifier_Predict']=np.argmax(global_pred,axis=1)

inp_df['Cumulative_Model_Score']=np.round((inp_df['Balanced_Classifier_Predict']+inp_df['Positively_Tuned_Classifier_Predict']+inp_df['Averaging_Classifier_Predict']+inp_df['Global_Classifier_Predict'])*100/4,2)
inp_df.drop(['Balanced_Classifier_Predict','Positively_Tuned_Classifier_Predict','Averaging_Classifier_Predict','Global_Classifier_Predict'],axis=1,inplace=True)

inp_df=inp_df[['temp_ind','Balanced_Classifier_Predict_Proba','Positively_Tuned_Classifier_Predict_Proba','Averaging_Classifier_Predict_Proba','Global_Classifier_Predict_Proba','Cumulative_Model_Score']]


final_df=pd.merge(init_df,inp_df,on=['temp_ind'],how='inner')
final_df.drop(['temp_ind'],axis=1,inplace=True)


final_df.to_csv(output_file,index=False,sep='\t')

