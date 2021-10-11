import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
import matplotlib.colors
from pathlib import Path
from Bio import SeqIO, SeqUtils
import gffpandas.gffpandas as gffpd
import seaborn as sns

from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.decomposition import PCA

def load_data(proj_dir):

    reference_df_path = proj_dir / 'Diel_LFC_results' / 'DGE_trained_reference_table.tsv'
    reference_df = pd.read_csv(reference_df_path, sep='\t', index_col='Gene ID')

    results_df_path = proj_dir / 'Diel_LFC_results' / 'DGE_trained_results_table.tsv'
    results_df = pd.read_csv(results_df_path, sep='\t', index_col='Gene ID')

    # ref_fasta = proj_dir / 'NATL2A_genome_references' / 'onlyNATL2A.fna'
    # ref_genome_record = SeqIO.read(ref_fasta, 'fasta')
    # ref_genome_seq = ref_genome_record.seq
    # gff_file = proj_dir / 'NATL2A_genome_references' / 'onlyNATL2A.gff'
    # annotation = gffpd.read_gff3(gff_file)
    # attributes_df = annotation.attributes_to_columns()

    return reference_df, results_df


def clustering_analysis(X, num_clusters):

    cluster_object = AgglomerativeClustering(n_clusters=num_clusters, compute_distances=True)

    clustering = cluster_object.fit(X)

    return clustering

def plot_diel_clusters(results_df, out_path, stacked=True):

    unique_clusters = results_df.index.get_level_values('cluster_label').unique()
    num_clusters = len(unique_clusters)

    if stacked:
        n_rows = 1
        n_cols = 1  
    else:
        n_cols = 5
        n_rows = num_clusters // n_cols
        if num_clusters % n_cols != 0: n_rows += 1

    cmap = plt.cm.viridis
    norm = matplotlib.colors.Normalize(vmin=min(unique_clusters), vmax=max(unique_clusters))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))
    plt.title('Layered Parental to Dark Adapted LFC Clusters ({} clusters)'.format(len(unique_clusters)))

    if stacked:
        ax = axes
        for cluster in unique_clusters:
            cluster_df = results_df.loc[cluster]
            cluster_arr = cluster_df.to_numpy()
            for r in cluster_arr:
                ax.plot(np.arange(len(r)), r, color=cmap(norm(cluster)), label=cluster, alpha=0.2)
        x_ticks = list(results_df.columns)
        ax.set_xticks(range(len(x_ticks)))
        ax.set_xticklabels(x_ticks)
    else:
        min_lfc = np.min(results_df.to_numpy())
        max_lfc = np.max(results_df.to_numpy())
        for i, ax in enumerate(axes.flat):
            if i < len(unique_clusters):
                cluster_df = results_df.loc[unique_clusters[i]]
                cluster_arr = cluster_df.to_numpy()
                for r in cluster_arr:
                    ax.plot(np.arange(len(r)), r, color=cmap(norm(unique_clusters[i])), label=unique_clusters[i], alpha=0.2)
                x_ticks = list(results_df.columns)
                ax.set_xticks(range(len(x_ticks)))
                ax.set_xticklabels(x_ticks)
                ax.set_ylim(min_lfc, max_lfc)
    
    plt.savefig(out_path)
    plt.close()

def plot_PCA_clusters(results_df, out_path):

    pca = PCA(n_components=2)
    transform = pca.fit_transform(results_df.to_numpy())

    clusters = results_df.index.get_level_values('cluster_label')
    unique_clusters = clusters.unique()

    color_dict = {}
    color=iter(plt.cm.rainbow(np.linspace(0,1,len(unique_clusters))))
    for i, t in enumerate(unique_clusters):
        color_dict[t] = next(color)
 
    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    ax = fig.gca()
    ax.scatter(transform[:,0], transform[:,1], c=list(map(color_dict.get, clusters)))

    patches = [mpatches.Patch(color=c, label=l) for l, c in color_dict.items()]
    ax.legend(handles=patches)

    plt.savefig(out_path)
    plt.close()

def plot_heatmap(results_df, out_path):
    figure = sns.clustermap(results_df)
    figure.savefig(out_path)
    plt.close()


def best_clusters_silhouette(X, out_path):

    cluster_options = range(2, 100)
    cluster_silhouette_scores = []

    for n in cluster_options:
        clustering = clustering_analysis(X, num_clusters=n)
        silhouette_avg = silhouette_score(X, clustering.labels_)
        cluster_silhouette_scores.append(silhouette_avg)
    
    fig = plt.figure(figsize=(10, 5), constrained_layout=True)
    ax = fig.gca()
    ax.plot(cluster_options, cluster_silhouette_scores)

    [print(x, '\t', score) for x, score in zip(cluster_options, cluster_silhouette_scores)]

    plt.savefig(out_path)
    plt.close()


def main(proj_dir):

    figure_dir = proj_dir / 'figures'
    figure_dir.mkdir(parents=True, exist_ok=True)

    diel_plotting = figure_dir / 'diel_plotting'
    diel_plotting.mkdir(parents=True, exist_ok=True)

    pca_plotting = figure_dir / 'pca_plotting'
    pca_plotting.mkdir(parents=True, exist_ok=True)

    heat_map = figure_dir / 'heat_map'
    heat_map.mkdir(parents=True, exist_ok=True)

    silhouette_dir = figure_dir / 'silhouette_plot'
    silhouette_dir.mkdir(parents=True, exist_ok=True)

    dataframe_dir = figure_dir / 'cluster_references'
    dataframe_dir.mkdir(parents=True, exist_ok=True)

    reference_df, results_df = load_data(proj_dir)

    X = results_df.to_numpy()

    best_clusters_silhouette(X, silhouette_dir / 'ideal_clustering.png')

    for num_clusters in [9, 18, 36, 53]:
        
        clustering = clustering_analysis(X, num_clusters=num_clusters)
        labels = clustering.labels_

        reference_df = reference_df.reset_index()
        results_df = results_df.reset_index()
        
        reference_df['cluster_label'] = labels
        results_df['cluster_label'] = labels

        reference_df = reference_df.set_index(['cluster_label', 'Gene ID'])
        results_df = results_df.set_index(['cluster_label', 'Gene ID'])

        plot_diel_clusters(results_df, diel_plotting / 'diel_plot_stacked_{}.png'.format(num_clusters), stacked=True)
        plot_diel_clusters(results_df, diel_plotting / 'diel_plot_separate_{}.png'.format(num_clusters), stacked=False)
        plot_PCA_clusters(results_df, pca_plotting / 'PCA_plot_{}.png'.format(num_clusters))
        plot_heatmap(results_df, heat_map / 'heat_map_{}.png'.format(num_clusters))
        reference_df.to_csv(dataframe_dir / '{}_clusters_reference_df.tsv'.format(num_clusters), sep='\t')
    

# # Parsons iMac path:
# proj_dir = Path(r'/Users/konnorvonemster/Dropbox (MIT)/Lab Member Files/Current Lab Members/Konnor von Emster/2021_trained_diel_LFC_motif_search')

# macbook pro path:
proj_dir = Path(r'/Users/kve/Dropbox (MIT)/Lab Member Files/Current Lab Members/Konnor von Emster/2021_trained_diel_LFC_motif_search')

main(proj_dir)