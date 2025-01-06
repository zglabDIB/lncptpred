import pandas as pd
import numpy as np
from tqdm import tqdm
import sys


def read_lncrna_seq(seq_file):
    lncrna_file=seq_file
    handle=open(lncrna_file,'r')
    for a in handle.readlines():
        lncrna=a

    handle.close()
    return lncrna

def select_interactive_segments(file_name):
    df=pd.read_csv(file_name,sep='\t')
    df=df[df.Cumulative_Model_Score>50]
    df.reset_index(drop=True,inplace=True)
    return df
    

def preprocess_interactive_segments(df):
    threshold=0.5
    col_list=['Balanced_Classifier_Predict_Proba','Positively_Tuned_Classifier_Predict_Proba','Averaging_Classifier_Predict_Proba', 'Global_Classifier_Predict_Proba']   
    df['Aggregate_Segments_Score'] = df[df[col_list] >= threshold].mean(axis=1)
    df.drop(['Balanced_Classifier_Predict_Proba','Positively_Tuned_Classifier_Predict_Proba','Averaging_Classifier_Predict_Proba', 
             'Global_Classifier_Predict_Proba','Cumulative_Model_Score'],axis=1,inplace=True)
    return df
    

def merge_overlapping_ranges(df):
    range_list=[]
    for i in tqdm(df.index.values):
        range_list.append((df['Window_Start'][i],df['Window_End'][i]))
    # Sort ranges based on the starting position
    sorted_ranges = sorted(range_list)

    merged_ranges = []
    current_start, current_end = sorted_ranges[0]

    for next_start, next_end in tqdm(sorted_ranges[1:]):
        if next_start <= current_end:  # Ranges overlap or are adjacent
            current_end = max(current_end, next_end)
        else:
            merged_ranges.append((current_start, current_end))
            current_start, current_end = next_start, next_end

    # Append the last range
    merged_ranges.append((current_start, current_end))

    return merged_ranges


def select_overlapping_segments(overlap_range,df):
    start_index=df[df.Window_Start==overlap_range[0]].index.values[0]
    end_index=df[df.Window_End==overlap_range[1]].index.values[0]
    df_new=df.iloc[start_index:end_index+1]
    return df_new



def final_process(final_df,df,merged_ranges,seq):
    count=0
    for overlap_range in tqdm(merged_ranges):
        df_new=select_overlapping_segments(overlap_range,df)
        final_df.loc[count,'Window_Start']=overlap_range[0]
        final_df.loc[count,'Window_End']=overlap_range[1]
        final_df.loc[count,'lncRNA_Sequence']=seq[overlap_range[0]-1:overlap_range[1]-1]
        final_df.loc[count,'Strand']=df.Strand.unique()[0]
        final_df.loc[count,'Protein_Name']=df.Protein_Name.unique()[0]
        final_df.loc[count,'Final_Interacting_Score']=df_new.Aggregate_Segments_Score.mean()
        count+=1
    final_df['Final_Interacting_Score']=final_df.Final_Interacting_Score.astype(np.float32).round(3)
    return final_df





seq_file=sys.argv[1]
intermediate_out_file=sys.argv[2]
output_file=sys.argv[3]

seq=read_lncrna_seq(seq_file)
final_df=pd.DataFrame(columns=['Window_Start', 'Window_End', 'lncRNA_Sequence', 'Strand',
       'Protein_Name','Final_Interacting_Score'])



df=select_interactive_segments(intermediate_out_file)
if(df.shape[0]>0):
    df=preprocess_interactive_segments(df)
    merged_ranges = merge_overlapping_ranges(df)
    final_df=final_process(final_df,df,merged_ranges,seq)
final_df.to_csv(output_file,index=False,sep='\t')


