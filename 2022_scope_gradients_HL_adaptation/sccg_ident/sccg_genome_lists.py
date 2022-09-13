import pandas as pd

df = pd.read_table("/nobackup1b/users/chisholmlab/img_proportal/data/lists/img_genome_lookup_table.txt")

print(df['type'].unique())
print(df['phylogeny'].unique())

df[(df['type']=='ISOLATE') & (df['phylogeny']=='Prochlorococcus')].to_csv("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/sccg_ident/pro_genomes_from_engaging.tsv", sep='\t')
df[(df['type']=='ISOLATE') & (df['phylogeny']=='Synechococcus')].to_csv("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/sccg_ident/syn_genomes_from_engaging.tsv", sep='\t')
