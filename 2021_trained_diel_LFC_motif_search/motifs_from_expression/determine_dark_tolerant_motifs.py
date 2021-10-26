import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
import matplotlib.colors
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffpandas.gffpandas as gffpd
import seaborn as sns
import time, sys
from gimmemotifs.denovo import gimme_motifs

import scipy.stats as stats

from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

TRANSLATION_TABLE = 11

def load_data(proj_dir):

    reference_df_path = proj_dir / 'Diel_LFCs' / 'DGE_trained_reference_table.tsv'
    reference_df = pd.read_csv(reference_df_path, sep='\t', index_col='Gene ID')

    results_df_path = proj_dir / 'Diel_LFCs' / 'DGE_trained_results_table.tsv'
    results_df = pd.read_csv(results_df_path, sep='\t', index_col='Gene ID')

    ref_fasta = proj_dir / 'NATL2A_genome_references' / 'onlyNATL2A.fna'
    ref_genome_record = SeqIO.read(ref_fasta, 'fasta')
    ref_genome_seq = ref_genome_record.seq

    gff_file = proj_dir / 'NATL2A_genome_references' / 'onlyNATL2A_w_sRNA.gff'
    annotation = gffpd.read_gff3(gff_file)
    attributes_df = annotation.attributes_to_columns()
    attributes_df['Gene ID'] = attributes_df['ID']
    attributes_df = attributes_df.set_index('Gene ID')

    reference_df = reference_df.join(attributes_df, on='Gene ID')

    return reference_df, results_df, ref_genome_seq


def clustering_analysis(X, num_clusters):

    cluster_object = AgglomerativeClustering(n_clusters=num_clusters, compute_distances=True)

    clustering = cluster_object.fit(X)

    return clustering

def plot_diel_clusters(results_df, out_path, stacked=True):

    unique_clusters = results_df.index.get_level_values('cluster_label').unique()
    num_clusters = len(unique_clusters) + 1

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
    fig.suptitle(f"LFC Between Parental and Dark Adapted Strains ({len(unique_clusters)} clusters)", size="xx-large")

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
                ax.set_title(f"Cluster {unique_clusters[i]}")
                cluster_df = results_df.loc[unique_clusters[i]]
                cluster_arr = cluster_df.to_numpy()
                columns = list(results_df.columns)

                x_points = np.arange(len(columns))

                for r in cluster_arr:
                    ax.plot(x_points, r, color=cmap(norm(unique_clusters[i])), alpha=0.2)

                sig_diffs = []
                sig_LFCs = []
                for col in columns:
                    lfc_arr = cluster_df[col].to_numpy()
                    avg = lfc_arr.mean()
                    std = lfc_arr.std()
                    z_val = abs(avg/std)
                    p_val = stats.norm.sf(z_val)*2
                    sig_diffs.append(p_val < 0.05)
                    sig_LFCs.append(abs(avg) > 1)

                y_range = max_lfc - min_lfc
                
                for x, sig_diff, sig_LFC in zip(x_points, sig_diffs, sig_LFCs):
                    if sig_LFC:
                        ax.plot(x, min_lfc + y_range*0.1, 'b*', label="Mean LFC > 1")
                    if sig_diff:
                        ax.plot(x, min_lfc + y_range*0.05, 'r*', label="DE cluster at 5% FDR")


                ax.set_xticks(x_points)
                ax.set_xticklabels(columns)
                ax.set_ylim(min_lfc, max_lfc)
            else:
                if i == len(unique_clusters):
                    legend_elements = [ Line2D([0], [0], marker='*', color='w', label="Mean LFC > 1", markerfacecolor='b', markersize=10),
                                        Line2D([0], [0], marker='*', color='w', label="DE cluster at 5% FDR", markerfacecolor='r', markersize=10)]
                    ax.legend(handles=legend_elements, loc=2)
                    ax.set_axis_off()
                else:
                    fig.delaxes(ax)
    
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
    for i, j in zip(cluster_options, cluster_silhouette_scores):
        print(i, '\t', j)

    plt.savefig(out_path)
    plt.close()


def main(proj_dir, num_clusters):

    results_dir = proj_dir / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)

    diel_plotting = results_dir / 'diel_plotting'
    diel_plotting.mkdir(parents=True, exist_ok=True)

    pca_plotting = results_dir / 'pca_plotting'
    pca_plotting.mkdir(parents=True, exist_ok=True)

    heat_map = results_dir / 'heat_map'
    heat_map.mkdir(parents=True, exist_ok=True)

    silhouette_dir = results_dir / 'silhouette_plot'
    silhouette_dir.mkdir(parents=True, exist_ok=True)

    dataframe_dir = results_dir / 'cluster_references'
    dataframe_dir.mkdir(parents=True, exist_ok=True)

    motif_fastas = results_dir / 'motif_fastas'
    motif_fastas.mkdir(parents=True, exist_ok=True)

    gimme_results = results_dir / 'gimme_motifs_results'
    gimme_results.mkdir(parents=True, exist_ok=True)

    reference_df, results_df, ref_genome_seq = load_data(proj_dir)

    X = results_df.to_numpy()

    # best_clusters_silhouette(X, silhouette_dir / 'ideal_clustering.png')

    # [9, 18, 39, 59]
    if num_clusters:
        
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

        # for cluster in set(labels):
        #     n_clust_dir = motif_fastas / '{:02}_clusters'.format(num_clusters)
        #     n_clust_dir.mkdir(parents=True, exist_ok=True)

        #     start_time = time.time()
        #     print(f"cluster {cluster:02} started of {num_clusters:02}")

        #     fasta_out = n_clust_dir / 'cluster_{:02}.fasta'.format(cluster)

        #     cluster_df = reference_df.loc[cluster]

        #     savePromoterFasta(ref_genome_seq, cluster_df, fasta_out)

        #     n_cluster_motif_dir = gimme_results / '{:02}_clusters'.format(num_clusters) / 'cluster_{:02}'.format(cluster)
        #     n_cluster_motif_dir.mkdir(parents=True, exist_ok=True)

        #     genome_path = proj_dir / 'NATL2A_genome_references' / 'onlyNATL2A.fna'

        #     try:
        #         gimme_motifs(str(fasta_out), str(n_cluster_motif_dir), params={"genome":str(genome_path)})
        #     except:
        #         pass
            
        #     td = time.time() - start_time
        #     print(f"cluster {cluster} finished of {num_clusters:02} in {td // 60} minutes")
            
    

def savePromoterFasta(ref_seq, gff_df, fasta_outpath):
    yfr_seq_records = []
    
    for gene in [gff_df.iloc[x] for x in range(len(gff_df))]:
        label = gene.ID
        start = gene.start
        end = gene.end
        strand = gene.strand
        product = gene['product']

        cds_sequence = ref_seq[start-1:end]
        if strand == '-':
            cds_sequence = cds_sequence.reverse_complement()

        # A good control I found was to translate the protein and observe if it made sense
        protein_sequence = cds_sequence.translate(table=TRANSLATION_TABLE)

        promoter_start = promoter_stop = 0

        if strand == '+':
            promoter_start = start-100
            promoter_stop = start+50
            promoter_region_sequence = ref_seq[promoter_start:promoter_stop]
        elif strand == '-':
            promoter_start = end-49
            promoter_stop = end+101
            promoter_region_sequence = ref_seq[promoter_start:promoter_stop]
            promoter_region_sequence = promoter_region_sequence.reverse_complement()

        # print(label, '\t', start, '\t', end, '\t', strand, '\t', cds_sequence[0:5] + '...' + cds_sequence[-5:])

        gene_seq_record = SeqRecord(promoter_region_sequence, 
                                    id=label, name=gene.name, 
                                    description='promoter region of {} making {} on strand {}'.format(label, product, strand))
        
        # print(gene_seq_record)
        yfr_seq_records.append(gene_seq_record)

    SeqIO.write(yfr_seq_records, fasta_outpath, "fasta")

# proj_dir = Path(sys.argv[1])

# num_clusters = int(sys.argv[2])

for num_clusters in [9, 18, 39, 59]:
    main(Path("/pool001/kve/2021_trained_diel_LFC_motif_search"), num_clusters)

