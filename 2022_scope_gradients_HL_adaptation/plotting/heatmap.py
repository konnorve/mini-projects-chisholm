from pathlib import Path
import pandas as pd
import numpy as np; np.set_printoptions(suppress=True)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns; sns.set()
from sklearn.decomposition import PCA

print('finished loading modules')

def plot_heatmap(results_df, cycog_df, out_path, cluster_cols=False, cluster_rows=False):

    stress_types = cycog_df["stress_type"].unique()
    stress_type_pal = sns.color_palette("tab10", len(stress_types))
    stress_type_color_dict = dict(zip(stress_types, stress_type_pal))
    stress_type_cols = results_df.columns.map(cycog_df['stress_type'].to_dict())
    stress_type_colors = stress_type_cols.map(stress_type_color_dict)

    stress_levels = ['Low', 'Medium', 'High']
    stress_level_pal = sns.light_palette("red", n_colors=6)
    stress_level_color_dict = dict(zip(stress_levels, stress_level_pal[1:4]))
    stress_level_cols = results_df.columns.map(cycog_df['stress_level'].to_dict())
    stress_level_colors = stress_level_cols.map(stress_level_color_dict)

    latitudes = results_df.index
    latitude_pal = sns.light_palette("blue", n_colors=len(latitudes), reverse=True)
    
    figure = sns.clustermap(results_df, figsize=(30,10), col_cluster=cluster_cols, row_cluster=cluster_rows, col_colors=[stress_type_colors, stress_level_colors], row_colors=[latitude_pal], dendrogram_ratio=0.15)

    stress_type_patches = [mpatches.Patch(facecolor=c, label=l) for l, c in stress_type_color_dict.items()]
    stress_level_patches = [mpatches.Patch(facecolor=c, label=l) for l, c in stress_level_color_dict.items()]
    
    legend_types = plt.legend(handles=stress_type_patches, title='Stress Types',
                    bbox_to_anchor=(0.9, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    
    legend_stress = plt.legend(handles=stress_level_patches, title='Stress Levels',
                    bbox_to_anchor=(0.95, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

    plt.gca().add_artist(legend_types)
    plt.gca().add_artist(legend_stress)

    figure.savefig(out_path)
    plt.close()

df = pd.read_table("scope-gradients-metag-highlightpro-stressgenes.tsv")
cycog_lookup_table = pd.read_table("cycog_lookup_table.tsv", index_col='gene_name')

print(cycog_lookup_table)

info_cols = ["Cruise", "library_id", "sample_id", "depth[m]", "latitude", "longitude", "time_local", "time_zulu"]

cycog_cols = [col for col in df.columns if 'CyCOG' in col]

def normalize(col):
    if col.max() > 1:
        return col / col.max()
    else:
        return col

for cruise in df["Cruise"].unique():
    print(cruise)
    cruise_df = df[df['Cruise']==cruise]
    cruise_subset = cruise_df[cycog_cols+["latitude"]]
    cruise_subset = cruise_subset.set_index("latitude")
    cruise_subset = cruise_subset.sort_index(ascending=False)
    cruise_subset = cruise_subset.apply(lambda col: normalize(col)).fillna(0)

    cruise_subset = cruise_subset.rename(columns={a:b for a,b in zip(cycog_lookup_table['cycog_id'], cycog_lookup_table.index)})

    print(cruise_subset)

    plot_heatmap(cruise_subset, cycog_lookup_table, f"{cruise} heatmap no clusters.png", cluster_cols=False, cluster_rows=False)
    plot_heatmap(cruise_subset, cycog_lookup_table, f"{cruise} heatmap row clusters.png", cluster_cols=False, cluster_rows=True)
    plot_heatmap(cruise_subset, cycog_lookup_table, f"{cruise} heatmap all clusters.png", cluster_cols=True, cluster_rows=True)


