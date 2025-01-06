import os

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

window_len_inp=input('Select Window length of lncRNA Sequences (Min:10 Max:30 Default:20)::')
if window_len_inp=='':
    window_len=20
else:
    window_len=int(window_len_inp)
    
if window_len<10:
    window_len=10
    
if window_len>30:
    window_len=30

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



    protein_list=["YBX1","ALKBH5","G3BP1","G3BP2","SF3B1","C17orf85","LIN28B","ZCCHC4","HuR","ZFP36","SRSF2","SRSF1","FMR1","SNRPA1","TAF15","EWSR1","ESR1","CAPRIN1","C22orf28","MSI2","LIN28A","RC3H1","YTHDC1","QKI","ZC3H7B","YTHDF1","YTHDF2","FUS","AUF1","METTL14","PTBP2","METTL3","WTAP"]
    print('Select Number below representing Protein Name: Default is HuR\n')
    print('''0=>YBX1\t\t1=>ALKBH5
            2=>G3BP1\t3=>G3BP2
            4=>SF3B1\t5=>C17orf85
            6=>LIN28B\t7=>ZCCHC4
            8=>HuR\t\t9=>ZFP36
            10=>SRSF2\t11=>SRSF1
            12=>FMR1\t13=>SNRPA1
            14=>TAF15\t15=>EWSR1
            16=>ESR1\t17=>CAPRIN1
            18=>C22orf28\t19=>MSI2
            20=>LIN28A\t21=>RC3H1
            22=>YTHDC1\t23=>QKI
            24=>ZC3H7B\t25=>YTHDF1
            26=>YTHDF2\t27=>FUS
            28=>AUF1\t29=>METTL14
            30=>PTBP2\t31=>METTL3
            32=>WTAP
    ''')
    protein_inp=input()

    if (protein_inp==''):
        protein_inp=-1
    else:
        protein_inp=int(protein_inp)

    if (protein_inp<0 or protein_inp>32):
        protein='HuR'
    else:
        protein=protein_list[protein_inp]


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
    os.remove(intermediate_out_file)