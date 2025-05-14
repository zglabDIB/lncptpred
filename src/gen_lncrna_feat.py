import pandas as pd 
import numpy as np
import itertools
from tqdm import tqdm
import re
import subprocess
from collections import Counter
import sys

input_file=sys.argv[1]
output_file=sys.argv[2]

# --- RNAfold runner ---
def run_rnafold(seq):
    result = subprocess.run(["RNAfold", "--noPS"], input=seq.encode(), capture_output=True, check=True)
    lines = result.stdout.decode().strip().split('\n')
    dotbracket = lines[1].split()[0]
    mfe = float(lines[1].split()[-1].strip('()'))
    return dotbracket, mfe

# --- Extract Secondary Structural Feature set ---
def extract_structure_features(dot_bracket):
    # Base pairs
    total_base_pairs = dot_bracket.count('(')

    # Unpaired bases
    unpaired_bases = dot_bracket.count('.')

    # Stems: count runs of base pairs like "(((" or ")))"
    stem_matches = re.findall(r'\(+', dot_bracket)
    stem_count = len(stem_matches)
    avg_stem_length = sum(len(match) for match in stem_matches) / stem_count if stem_count else 0.0

    # Hairpins: pattern (....)
    hairpins = re.findall(r'\([^()\.]{0,1}\.*[^()\.]{0,1}\)', dot_bracket)
    hairpin_count = len(hairpins)
    avg_hairpin_size = (
        sum(m.count('.') for m in hairpins) / hairpin_count if hairpin_count else 0.0
    )

    # Bulges: look for unpaired dots adjacent to a single bracket
    bulges = re.findall(r'\(\.+\(', dot_bracket) + re.findall(r'\)\.+\)', dot_bracket)
    bulge_count = len(bulges)

    # Multiloops: roughly defined as regions where 3 or more stems meet
    # A simple approximation: count regions with ")(." more than once
    multiloop_count = len(re.findall(r'\)[.]+\(', dot_bracket)) - 1
    multiloop_count = max(multiloop_count, 0)


    base_pair_density=round(dot_bracket.count('(') / len(dot_bracket), 4)

    return {
        'total_base_pairs': total_base_pairs,
        'unpaired_bases': unpaired_bases,
        'stem_count': stem_count,
        'avg_stem_length': round(avg_stem_length, 2),
        'hairpin_count': hairpin_count,
        'avg_hairpin_size': round(avg_hairpin_size, 2),
        'bulge_count': bulge_count,
        'multiloop_count': multiloop_count,
        'base_pair_density':base_pair_density
    }


# --- Generate Nucleotide Features ---
def gen_nucleotide_feat(df):
    sequences = df['lncRNA_Sequence'].tolist()
    
    def get_kmer_fractions(seq, k):
        total = len(seq) - k + 1
        if total <= 0:
            return {}
        kmers = [seq[i:i+k] for i in range(total)]
        counts = Counter(kmers)
        return {kmer: counts[kmer]/total for kmer in all_kmers[k]}
    
    # Precompute k-mers for k = 1 to 4
    all_kmers = {
        k: [''.join(p) for p in itertools.product('ACGT', repeat=k)]
        for k in range(1, 5)
    }
    
    features = []
    for seq in tqdm(sequences):
        feat = {}
        for k in range(1, 5):
            kmer_fractions = get_kmer_fractions(seq, k)
            for kmer in all_kmers[k]:
                feat[f'{kmer}_lncrna_fraction'] = kmer_fractions.get(kmer, 0.0)
        features.append(feat)
    
    feat_df = pd.DataFrame(features)
    result_df = pd.concat([df.reset_index(drop=True), feat_df], axis=1)
    
    return result_df



df=pd.read_csv(input_file,sep='\t')
df.loc[df['Strand']=='-','Strand']=-1
df.loc[df['Strand']=='+','Strand']=1
df['Strand']=df['Strand'].astype(np.int8)


features = []
for seq in tqdm(df["lncRNA_Sequence"]):
    dotbracket, mfe = run_rnafold(seq)
    feats = extract_structure_features(dotbracket)
    feats["mfe"] = mfe
    features.append(feats)

# --- Merge with original dataframe ---
features_df = pd.DataFrame(features)
df = pd.concat([df, features_df], axis=1)


df=gen_nucleotide_feat(df)
df.drop(['Window_Start','Window_End','lncRNA_Sequence'],inplace=True,axis=1)


cols_arrange=list(df.columns[1:])+[df.columns[0]]
df=df[cols_arrange]


df.to_csv(output_file,index=False,sep='\t')



