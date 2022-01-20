from pathlib import Path
import pandas as pd
import numpy as np
import re
import logging as log
import urllib.error

import dash
from dash import dcc
from dash import html
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output

import Bio.KEGG.REST as kegg_rest

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)


# columns of importance: 'BRITE', 'GOs', 'KEGG_Module', 'KEGG_Pathway', 'KEGG_Reaction', 'KEGG_rclass'

concat_df = pd.read_csv("concat_df_w_terms.tsv", sep='\t')
ref_table = pd.read_csv("kegg_ids.tsv", sep='\t')
    
concat_df_kegg_cols = ref_table.ref_col.unique()

concat_df = concat_df.fillna(value={col:"" for col in concat_df_kegg_cols})

app = dash.Dash(__name__)

@app.callback(
    Output('pathway-id', 'options'),
    Input('pathway-category', 'value'))
def set_pathway_options(pathway_category):
    ref_table_subset = ref_table[ref_table['ref_col']==pathway_category]
    ref_table_subset = ref_table_subset.dropna()
    ref_table_subset = ref_table_subset.sort_values('kegg_term')
    ref_options = [{'label':row.loc['kegg_term'], 'value':row.loc['kegg_id']} for i, row in ref_table_subset.iterrows()]
    return ref_options

@app.callback(
    Output('scatter-plot', 'figure'),
    Input('pathway-category', 'value'),
    Input('pathway-id', 'value'))
def update_scatter(pathway_category, pathway_id):
    fig = px.scatter(concat_df, 
                    x='log2FoldChange_exp11', 
                    y='log2FoldChange_exp4', 
                    hover_data=['Description', 'EC'], 
                    height = 1600, width=1600,
                    labels={
                        "log2FoldChange_exp11": "Effect of ProMS (LFC)",
                        "log2FoldChange_exp4": "Effect of co-culutre (LFC)"
                    },
                    title='Effects of ProMS and SAR11 on MIT9301 gene expression',
                    color=concat_df[pathway_category].str.contains(pathway_id)
                    )
                    

    fig.update_traces(marker = {
                            # 'color' : 'red',
                            'opacity' : 0.5,
                        },
                        selector=dict(mode='markers'))
    return fig


app.layout = html.Div([
    dcc.Dropdown(
        id='pathway-category',
        options=[{'label':val, 'value': val} for val in concat_df_kegg_cols],
        value='KEGG_Pathway'
    ),
    dcc.Dropdown(
        id='pathway-id'
        # value='map00195'
    ),
    dcc.Graph(id='scatter-plot')
])

if __name__ == "__main__":
    app.run_server(debug=True)