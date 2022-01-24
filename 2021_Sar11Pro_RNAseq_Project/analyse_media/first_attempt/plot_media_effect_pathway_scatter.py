from collections import namedtuple
from pathlib import Path
from matplotlib.colors import same_color
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

import re
import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.INFO)

def plot_LFCs(concat_df, include_control=True, cumulative=True, fig_out_path=None, alpha=1):

    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    ax = fig.gca()

    joined_df = concat_df.xs(key='log2FoldChange', axis=1, level=1)

    if include_control:
        joined_df['Control'] = 0

    for row in joined_df.to_numpy():
        if cumulative:
            ax.plot(range(len(row)), np.cumsum(row), alpha=alpha)
        else:
            ax.plot(range(len(row)), row, alpha=alpha)

    ax.set_xticks(range(len(row)))
    ax.set_xticklabels(joined_df.columns)
    
    ax.set_xlabel('addition')
    ax.set_ylabel('log2FoldChange')
    
    if fig_out_path:
        plt.savefig(fig_out_path)
    else:
        plt.show()
    plt.close()

def plot_LFCs_scatter(concat_df, fig_out_path=None, alpha=1):

    joined_df = concat_df.xs(key='log2FoldChange', axis=1, level=1)

    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    ax = fig.gca()

    for row in joined_df.to_numpy():
        ax.plot(row[0], row[1], "r.", alpha=alpha)
    
    if fig_out_path:
        plt.savefig(fig_out_path)
    else:
        plt.show()
    plt.close()

def plot_pathways_LFC_scatter(concat_df, bool_series, title, fig_out_path=None, alpha=1):

    joined_df = concat_df.xs(key='log2FoldChange', axis=1, level=1)

    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    ax = fig.gca()

    non_pathway_df = joined_df[~bool_series]
    for row in non_pathway_df.to_numpy():
        ax.plot(row[0], row[1], marker='.', color='gray', alpha=alpha)

    pathway_df = joined_df[bool_series]
    for row in pathway_df.to_numpy():
        ax.plot(row[0], row[1], marker='.', color='red', alpha=alpha)

    ax.set_title(title)
    ax.set_xlabel("LFC due to ProMS treatment")
    ax.set_ylabel("LFC due to Sar11 addition")
    
    if fig_out_path:
        plt.savefig(fig_out_path)
    else:
        plt.show()
    plt.close()

def plot_all_pathways(complex_df, go_dict, category_dict, out_dir):
    for col, readable_dict in go_dict.items():
        category = category_dict[col]
        category_dir = out_dir / category
        category_dir.mkdir(exist_ok=True)
        log.info(f"category: {category}")
        len_category = len(readable_dict)
        for i, (h, name) in enumerate(readable_dict.items()):
            out_path = category_dir /  f"{make_file_name(name)}.png"
            try:
                bool_series = complex_df['pathway', col + ' hash'].notna() & complex_df['pathway', col +  ' hash'].str.contains(h)
                if sum(bool_series) == 0:
                    log.error(f"No terms for {h} with {name}")
                log.info(f"{i:03}/{len_category:03} ({(i / len_category):2.2%}) \t {sum(bool_series)} \t {name}")
                plot_pathways_LFC_scatter(complex_df, bool_series, name, fig_out_path=out_path, alpha=0.5)
            except:
                log.error(f"Could not subset using hash {h} with {name}")
                raise


experiment_4_results_df = pd.read_csv("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/results/experiments/experiment_4/DGE_tables/experiment_4_MIT9301_DGE_all.tsv", index_col=0, sep="\t")
experiment_11_results_df = pd.read_csv("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/results/experiments/experiment_11/DGE_tables/experiment_11_MIT9301_DGE_all.tsv", index_col=0,  sep="\t")

experiment_4_results_df = experiment_4_results_df[~experiment_4_results_df.index.duplicated()]
experiment_11_results_df = experiment_11_results_df[~experiment_11_results_df.index.duplicated()]

concat_df = pd.concat({"proMS" : experiment_11_results_df, "co-cultures" : experiment_4_results_df}, axis=1, names=['treatment'])

pathway_df = pd.read_csv("/home/kve/scripts/mini_projects/2021_Sar11Pro_RNAseq_Project/analyse_media/pathway_tools_df_MIT9301.tsv", sep='\t')
pathway_df = pathway_df.set_index("NCBI-Protein")
pathway_df = pd.concat({'pathway':pathway_df}, axis=1)
pathway_df = pathway_df[pathway_df.index.notna() & ~pathway_df.index.duplicated()]

complex_df = concat_df.join(pathway_df, on=[('proMS','protein_id')])




def make_file_name(s):
    s = re.sub(r'<.*?>', '', s).replace("/", " or ")
    s = re.sub('[\&\';.\[\]+]', '', s).replace(" ", "_")
    return s

def make_readable(s):
    s = re.sub(r'<.*?>', '', s)
    s = re.sub('[&;.]', '', s)
    return s

go_dict = {}
for col in ['Pathways of enzyme', 'GO terms (biological process)', 'GO terms (cellular component)', 'GO terms (molecular function)']:
    gene_group_terms = complex_df.loc[complex_df['pathway',col].notna()]['pathway', col].to_list()
    gene_group_terms_list = [s.split(" // ") for s in gene_group_terms]
    gene_group_terms_set = {item for sublist in gene_group_terms_list for item in sublist}

    hash_dict = {}
    readable_dict = {}
    for term in gene_group_terms_set:
        h = hex(hash(term))
        hash_dict[term] = h
        readable_dict[h] = make_readable(term)

    gene_group_terms_hash_list = [[hash_dict[term] for term in gene_terms_list] for gene_terms_list in gene_group_terms_list]
    gene_hash_list = [" // ".join(l) for l in gene_group_terms_hash_list]

    gene_group_terms_readable_list = [[readable_dict[hash_dict[term]] for term in gene_terms_list] for gene_terms_list in gene_group_terms_list]
    gene_readable_list = [" // ".join(l) for l in gene_group_terms_readable_list]

    complex_df.loc[complex_df['pathway',col].notna(), ('pathway', col + " hash")] = gene_hash_list
    complex_df.loc[complex_df['pathway',col].notna(), ('pathway', col + " readable")] = gene_readable_list

    go_dict[col] = readable_dict

scatter_plot_path = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/analyse_media/pathway_scatter_plots")

complex_df.to_csv(scatter_plot_path / "media_co-culture_treatments_w_pathways.tsv", sep='\t')

category_dict = {col:re.search(r'\((.*?)\)', col).group(1).replace(" ", "_") for col in ['GO terms (biological process)', 'GO terms (cellular component)', 'GO terms (molecular function)']}
category_dict['Pathways of enzyme'] = 'pathwaytools_enzyme_pathway'

plot_all_pathways(complex_df, go_dict, category_dict, scatter_plot_path)