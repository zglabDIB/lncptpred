import os
import pandas as pd

lncrna_seq_file='../input/'+'lncrna_seq.txt';
initial_file='../input/'+'lncrna_prot_divide_seq.txt';
lncrna_file='../input/'+'lncrna_gen_feat.txt';
prot_file='../input/'+'prot_gen_feat.txt';
feat_file='../input/'+'final_feature.txt';
intermediate_out_file='../output/'+'intermediate_output.txt';
out_file='../output/'+'prediction_output.txt';
img_file='../output/'+'motif_logo.png';


handle=open(lncrna_seq_file,'r')
for a in handle.readlines():
    lncrna=a

handle.close()

total_len=len(lncrna)

window_len_inp=input('Select Window length of lncRNA Sequences (Min:10 Max:40 Default:20)::')
if window_len_inp=='':
    window_len=20
else:
    window_len=int(window_len_inp)
    
if window_len<10:
    window_len=10
    
if window_len>40:
    window_len=40

if (total_len<window_len):
    print('Input Sequence Length must be greater or equal to window length!!!!!!!!!!!!!!!!!!!')

else:
    shift_size_inp=input('Select Shifting size of lncRNA Sequence Windows (Min:1 Max:5 Default:1)::')
    if shift_size_inp=='':
        shift_size=1
    else:
        shift_size=int(shift_size_inp)
        
    if shift_size<1:
        shift_size=1
        
    if shift_size>5:
        shift_size=5



    strand_inp=input('Select Strands (+ or -   Default +)::')
    if strand_inp=='-':
        strand='-'
    else:
        strand='+'



    print('Select Number below representing Protein Name: Default is HuR\n')
    prot_df=pd.read_csv('../dataset/Final_Protein.txt',sep='\t')
    prot_dict = dict(zip(prot_df['Prot_ID'], prot_df['Protein']))

    count = 0
    for k, v in prot_dict.items():
        print(f"{k}==>{v}", end='\t')
        count += 1
        if count % 8 == 0:
            print()  # Newline after every 4 pairs



    protein_inp=input()

    if (protein_inp==''):
        protein_inp=-1
    else:
        protein_inp=int(protein_inp)

    if (protein_inp<1 or protein_inp>88):
        protein='HuR'
    else:
        protein=prot_dict[protein_inp]


    print('Protein Name:'+protein)
    print('Window Length:'+str(window_len))
    print('Shift Size:'+str(shift_size))
    print('Strand:'+strand)


    command_1="python split_seq.py "+lncrna_seq_file+" "+str(window_len)+" "+str(shift_size)+" "+strand+" "+protein+" "+initial_file
    command_2="python gen_lncrna_feat.py "+initial_file+" "+lncrna_file
    command_3="python gen_protein_feat.py "+protein+" "+prot_file
    command_4="python merged_features.py "+lncrna_file+" "+prot_file+" "+feat_file
    command_5="python intermediate_interactions.py "+initial_file+" "+feat_file+" "+intermediate_out_file
    command_6="python predict_interactions.py "+lncrna_seq_file+" "+intermediate_out_file+" "+out_file
    command_7="python gen_logo.py "+intermediate_out_file+" "+img_file+" "+str(window_len)



    os.system(command_1)
    os.system(command_2)
    os.system(command_3)
    os.system(command_4)
    os.system(command_5)
    os.system(command_6)
    os.system(command_7)

    os.remove(initial_file)
    os.remove(lncrna_file)
    os.remove(prot_file)
    os.remove(feat_file)
