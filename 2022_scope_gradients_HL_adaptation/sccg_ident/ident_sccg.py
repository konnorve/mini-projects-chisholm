import pandas as pd
from pathlib import Path
import plotly.express as px

wd = Path('/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/sccg_ident')
pro_genomes = pd.read_table(wd / 'pro_sccg_genomes_from_engaging.tsv')
syn_genomes = pd.read_table(wd / 'syn_sccg_genomes_from_engaging.tsv')

cycog_genes = pd.read_table('/nobackup1b/users/chisholmlab/img_proportal/cycogs/CyCOGs-v6.0/cycog-genes.tsv', index_col='gene_iid')

assert pro_genomes.columns.to_list() == syn_genomes.columns.to_list()

print(f"genome columns: {pro_genomes.columns}")
print(f"cycog_genes columns: {cycog_genes.columns}")


pro_genomes_oid = set(pro_genomes.img_genome_id)
syn_genomes_oid = set(syn_genomes.img_genome_id)
cycog_oids = set(cycog_genes.taxon_oid)
pro_oid_intersection = pro_genomes_oid.intersection(cycog_oids)
pro_oid_difference = pro_genomes_oid.difference(cycog_oids)
syn_oid_intersection = syn_genomes_oid.intersection(cycog_oids)
syn_oid_difference = syn_genomes_oid.difference(cycog_oids)

for s in [pro_oid_intersection, pro_oid_difference, syn_oid_intersection, syn_oid_difference]:
    print(len(s), s)

for s in [pro_genomes, syn_genomes, cycog_genes]:
    print(len(s))

pro_genomes = pro_genomes[pro_genomes['img_genome_id'].isin(pro_oid_intersection)]
syn_genomes = syn_genomes[syn_genomes['img_genome_id'].isin(syn_oid_intersection)]
pro_cycogs = cycog_genes[cycog_genes['taxon_oid'].isin(pro_oid_intersection)]
syn_cycogs = cycog_genes[cycog_genes['taxon_oid'].isin(syn_oid_intersection)]

for s in [pro_genomes, syn_genomes, pro_cycogs, syn_cycogs]:
    print(len(s))

# cycogs as rows, genomes as columns
#subset cycogs on 
pro_cycog_occurrences = pd.get_dummies(pro_cycogs[['taxon_oid', 'cycog_iid']], columns=['taxon_oid']).groupby(['cycog_iid']).sum()
pro_cycog_occurrences = pro_cycog_occurrences[pro_cycog_occurrences <= 1].dropna()
pro_cycog_pct = pro_cycog_occurrences.sum(axis=1) / len(pro_cycog_occurrences.columns)
pro_cycog_pct_value_counts = pro_cycog_pct.value_counts().sort_index()

# fig = px.line(x=pro_cycog_pct_value_counts.cumsum().values, y=pro_cycog_pct_value_counts.index)
# fig.write_html('pro_cycog_pct_value_counts.html')
# fig.write_image('pro_cycog_pct_value_counts.png')


syn_cycog_occurrences = pd.get_dummies(syn_cycogs[['taxon_oid', 'cycog_iid']], columns=['taxon_oid']).groupby(['cycog_iid']).sum()
syn_cycog_occurrences = syn_cycog_occurrences[syn_cycog_occurrences <= 1].dropna()
syn_cycog_pct = syn_cycog_occurrences.sum(axis=1) / len(syn_cycog_occurrences.columns)
syn_cycog_pct_value_counts = syn_cycog_pct.value_counts().sort_index()

# fig = px.line(x=syn_cycog_pct_value_counts.cumsum().values, y=syn_cycog_pct_value_counts.index)
# fig.write_html('syn_cycog_pct_value_counts.html')
# fig.write_image('syn_cycog_pct_value_counts.png')


pro_sccg_cycogs = pro_cycog_pct[pro_cycog_pct > 0.95].index.to_list()
syn_sccg_cycogs = syn_cycog_pct[syn_cycog_pct > 0.95].index.to_list()

[print(len(l)) for l in [pro_sccg_cycogs, syn_sccg_cycogs]]

# with open('pro_sccg.txt', 'w') as f:
#     f.write("\n".join(pro_sccg_cycogs))

# with open('syn_sccg.txt', 'w') as f:
#     f.write("\n".join(syn_sccg_cycogs))
