import pandas as pd
from pathlib import Path


all_path = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/3_flagstat/G1St3100m_S30.tsv")
paired_path = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/3_flagstat/G1St3100m_S30_paired.tsv")
unpaired_path = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/3_flagstat/G1St3100m_S30_unpaired.tsv")

all_results = pd.read_table(all_path, names=['all', 'qc-fail', 'ident'], usecols=['all', 'ident']).set_index('ident')
paired_results = pd.read_table(paired_path, names=['paired', 'qc-fail', 'ident'], usecols=['paired', 'ident']).set_index('ident')
unpaired_results = pd.read_table(unpaired_path, names=['unpaired', 'qc-fail', 'ident'], usecols=['unpaired', 'ident']).set_index('ident')

df = pd.concat([all_results, paired_results, unpaired_results], axis=1)

subsetter = pd.concat([df[col].str.isnumeric() for col in df.columns], axis=1).all(axis=1)

df = df[subsetter].astype('int32')

print(df)

df['sum'] = df['paired'] + df['unpaired']

df['equal'] = df['sum'] == df['all']

print(df)

