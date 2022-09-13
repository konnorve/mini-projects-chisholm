import pandas as pd
from pathlib import Path

wd = Path('/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/sccg_ident/')

dfs = {}
for table in wd.glob("*.tsv"):
    dfs[table.stem] = pd.read_table(table)

# print("\n".join([f"{k}: {v.columns}" for k,v in dfs.items()]))

img_pro_genome_set_engaging = set(dfs['pro_genomes_from_engaging']['img_genome_id'].values)
img_syn_genome_set_engaging = set(dfs['syn_genomes_from_engaging']['img_genome_id'].values)
img_pro_genome_set_proportal = set(dfs['pro_sccg_genomes_from_IMG']['IMG Genome ID '].values)
img_syn_genome_set_proportal = set(dfs['syn_sccg_genomes_from_IMG']['IMG Genome ID '].values)

pro_intersection = img_pro_genome_set_engaging.intersection(img_pro_genome_set_proportal)
pro_engaging_only = img_pro_genome_set_engaging.difference(img_pro_genome_set_proportal)
pro_proportal_only = img_pro_genome_set_proportal.difference(img_pro_genome_set_engaging)

syn_intersection = img_syn_genome_set_engaging.intersection(img_syn_genome_set_proportal)
syn_engaging_only = img_syn_genome_set_engaging.difference(img_syn_genome_set_proportal)
syn_proportal_only = img_syn_genome_set_proportal.difference(img_syn_genome_set_engaging)

print(pro_intersection)
print()
print(pro_engaging_only)
print()
print(pro_proportal_only)
print()
print(syn_intersection)
print()
print(syn_engaging_only)
print()
print(syn_proportal_only)
print()