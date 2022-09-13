import pandas as pd
import yaml
import numpy as np
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

occurance_table_path = Path('/nfs/chisholmlab001/kve/2022_SNPs_Dark_Adapted_Genomes/analysis/variant_occurance_table.bcftools_standard.tsv')
treatment_yaml_path = Path('/nfs/chisholmlab001/kve/2022_SNPs_Dark_Adapted_Genomes/analysis/dark_adapted.yaml')

occurance_table = pd.read_table(occurance_table_path, index_col='id')
with open(treatment_yaml_path, 'r') as f:
    treatment_dict = yaml.safe_load(f)
treatment_df = pd.Series(treatment_dict)

control_occurances = occurance_table[treatment_df[treatment_df=='control'].index]
pheno_occurances = occurance_table[treatment_df[treatment_df=='pheno'].index]

control_rates = control_occurances.sum(axis=1) / len(control_occurances.columns)
pheno_rates = pheno_occurances.sum(axis=1) / len(pheno_occurances.columns)

control_rates.name = 'control_rates'
pheno_rates.name = 'pheno_rates'

df = pd.DataFrame([control_rates, pheno_rates]).transpose()

df['diff'] = df['control_rates'] - df['pheno_rates']

df = df.sort_values('diff').reset_index()

snp_occurances = pd.get_dummies(df[['control_rates', 'pheno_rates']], columns=['control_rates'], prefix=None).groupby(['pheno_rates']).sum()

fig = go.Figure()
fig.add_trace(
    go.Heatmap(
        z=np.log1p(snp_occurances.to_numpy()+1), 
        zsmooth = 'best',
        x=snp_occurances.columns,
        y=snp_occurances.index,
    )
)
fig.update_layout(
    xaxis_title='control_rates',
    yaxis_title='pheno_rates'
)
fig.write_html("control_vs_pheno_snp_rate_heatmap.html")
fig = px.scatter(df, x='control_rates', y='pheno_rates', hover_data=['id'])
fig.write_html("control_vs_pheno_snp_rate_scatter.html")