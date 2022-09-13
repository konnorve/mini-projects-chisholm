import pandas as pd
from pathlib import Path
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

img_genomes = pd.read_table("/nobackup1b/users/chisholmlab/img_proportal/data/lists/img_genome_lookup_table.txt")
cycog_genes = pd.read_table('/nobackup1b/users/chisholmlab/img_proportal/cycogs/CyCOGs-v6.0/cycog-genes.tsv', index_col='gene_iid')

# print(img_genomes)

# print(cycog_genes)

pro_sccg = set(Path('pro_sccg.txt').read_text().splitlines())
syn_sccg = set(Path('syn_sccg.txt').read_text().splitlines())

pro_sags = img_genomes[(img_genomes['type']=='SAG') & (img_genomes['phylogeny']=='Prochlorococcus')]
pro_sag_genome_ids = set(pro_sags.img_genome_id.to_list())
pro_sag_cycogs = cycog_genes[cycog_genes['taxon_oid'].isin(pro_sag_genome_ids)]
pro_sag_sccg_cycogs = pro_sag_cycogs[pro_sag_cycogs['cycog_iid'].isin(pro_sccg)]
pro_cycog_sag_pcts = pro_sag_sccg_cycogs.groupby('cycog_iid')['taxon_oid'].nunique() / pro_sag_cycogs['taxon_oid'].nunique()

fig = px.histogram(pro_cycog_sag_pcts.rename("proportion of sags carrying each putative sccg"))
mean = np.mean(pro_cycog_sag_pcts.values)
std = np.std(pro_cycog_sag_pcts.values)
for n in range(-3,4):
    x=mean+(n*std)
    fig.add_trace(go.Scatter(x=[x,x], 
                            y=[0,150], 
                            mode='lines', 
                            line_color='blue',
                            opacity=pow(2, -abs(n)),
                            name=f'mean {n:+} sd'))
fig.add_trace(go.Scatter(x=[0.606,0.606], y=[0,150], mode='lines', line_color='red', name='checkm avg completeness'))
fig.write_image("pro_cycog_sag_pcts.png", height=500, width=1200)

# print(pro_sag_cycogs)

syn_sags = img_genomes[(img_genomes['type']=='SAG') & (img_genomes['phylogeny']=='Synechococcus')]
syn_sag_genome_ids = set(syn_sags.img_genome_id.to_list())
syn_sag_cycogs = cycog_genes[cycog_genes['taxon_oid'].isin(syn_sag_genome_ids)]
syn_sag_sccg_cycogs = syn_sag_cycogs[syn_sag_cycogs['cycog_iid'].isin(syn_sccg)]
syn_cycog_sag_pcts = syn_sag_sccg_cycogs.groupby('cycog_iid')['taxon_oid'].nunique() / syn_sag_cycogs['taxon_oid'].nunique()

fig = px.histogram(syn_cycog_sag_pcts.rename("proportion of sags carrying each putative sccg"))
mean = np.mean(syn_cycog_sag_pcts.values)
std = np.std(syn_cycog_sag_pcts.values)
mean = np.mean(syn_cycog_sag_pcts.values)
std = np.std(syn_cycog_sag_pcts.values)
for n in range(-3,4):
    x=mean+(n*std)
    fig.add_trace(go.Scatter(x=[x,x], 
                            y=[0,150], 
                            mode='lines', 
                            line_color='blue',
                            opacity=pow(2, -abs(n)),
                            name=f'mean {n:+} sd'))
fig.add_trace(go.Scatter(x=[0.48,0.48], y=[0,150], mode='lines', line_color='red', name='checkm avg completeness'))
fig.write_image("syn_cycog_sag_pcts.png", height=500, width=1200)

# print(syn_sag_cycogs)

# in_sccg_not_in_sags = pro_sccg.difference(set(pro_sag_cycogs.cycog_iid.to_list()))
# print(in_sccg_not_in_sags)

# in_sccg_not_in_sags = syn_sccg.difference(set(syn_sag_cycogs.cycog_iid.to_list()))
# print(in_sccg_not_in_sags)