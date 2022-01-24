
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import gffpandas.gffpandas as gffpd
import unicodedata
import re

results_dir = Path("/nfs/chisholmlab001/kve/2021_dark_adapted_transcriptome/results/")

experiment_dir = results_dir / "experiments"
plotting_dir = results_dir / 'gene_expression_plots'

log_count_table = experiment_dir / "experiment_all" / "DEseq2out" / "experiment_all_NATL2A_rlog.tsv"
rlog_df = pd.read_csv(log_count_table, sep='\t', index_col="long_ID")

rlog_df.columns = pd.MultiIndex.from_product([['Control', 'Pheno'], [0, 4, 8, 13, 16, 20, 24], [1, 2, 3]], names=['treatment', 'time', 'replicate'])

rlog_df = rlog_df[~rlog_df.index.duplicated()]

rlog_mean_df = rlog_df.groupby(level=['treatment', 'time'], axis='columns').mean()

rlog_std_df = rlog_df.groupby(level=['treatment', 'time'], axis='columns').std()


dfs = []
for time in [0, 4, 8, 13, 16, 20, 24]:
    result_path = experiment_dir / f"experiment_{time}" / "DGE_tables" / f"experiment_{time}_NATL2A_DGE_all.tsv"
    result_df = pd.read_csv(result_path, sep="\t", index_col="long_ID")

    # annotation_data = result_df[["product", "protein_id", "locus_tag"]]
    # annotation_data = annotation_data.drop_duplicates()

    significance_data = result_df[["log2FoldChange", "padj"]]
    significance_data.columns = pd.MultiIndex.from_product([[time], significance_data.columns], names=['time', 'significance'])

    significance_data = significance_data.drop_duplicates()

    dfs.append(significance_data)

result_df = pd.concat(dfs, axis=1)

result_df = result_df.swaplevel("time", "significance", axis=1)

# conversion_df = pd.read_csv("gene_label_conversion_table.tsv", sep='\t', index_col='gene')

# annotation_data = annotation_data.join(conversion_df, how='left')


gff_path = Path("/nfs/chisholmlab001/kve/2021_dark_adapted_transcriptome/input_data/culture_genome_annotations/NATL2A.gff")

annotation = gffpd.read_gff3(gff_path)
annotation_data = annotation.attributes_to_columns()

conversion_df = pd.read_csv("natl2a_convertion_table2.tsv", sep='\t')
conversion_df["locus_tag"] = conversion_df['NCBI ID_2']
conversion_df = conversion_df[conversion_df["locus_tag"].notna()]

annotation_data = annotation_data.merge(conversion_df, how='left', on="locus_tag")
annotation_data = annotation_data.set_index("ID")

annotation_data = annotation_data[annotation_data['type'].isin(['sRNA', 'CDS'])]
annotation_data = annotation_data.drop_duplicates()

annotation_data

# annotation_data["Gene Name"].to_list()
clock_renaming =   {'cds-WP_011294127.1': 'RpaA',
                    'cds-WP_011294763.1': 'SasA',
                    'cds-WP_011295261.1': 'KaiC',
                    'cds-WP_011295262.1': 'KaiB',
                    'cds-WP_011295475.1': 'LdpA'}
                    
for k, v in clock_renaming.items():
    annotation_data.loc[k, "Genbank Annotation"] = f"{v} gene"

def slugify(value, allow_unicode=False):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')


@mpl.rc_context({
    'lines.linewidth': 6, 
    'lines.marker':'o', 
    'lines.markersize':18, 
    'legend.fontsize': 'x-large',
    'axes.labelsize': 'x-large',
    'axes.titlesize':'x-large',
    'xtick.labelsize':'x-large',
    'ytick.labelsize':'x-large'})
def plot_gene(ax, gene_ID, rlog_mean_df, rlog_std_df, results_df, annotation_data, night_periods, night_color, attr_dict):

    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    mean_series = rlog_mean_df.loc[gene_ID]
    std_series = rlog_std_df.loc[gene_ID]

    for treatment in ['Control', 'Pheno']:
        for (s, e) in night_periods:
            ax.axvspan(s, e, color=night_color)
        ax.errorbar(x=mean_series[treatment].index, y=mean_series[treatment], yerr=std_series[treatment], capsize=6, capthick=3, label=treatment, color=attr_dict[treatment]['color'])
    
    # axes limits
    bottom, top = ax.get_ylim()
    y_range = top - bottom
    ax.set_ylim(top-(y_range*1.1), top)
    ax.set_xlim(-1, 25)
    
    #plotting significance
    try:
        significance_df = result_df.loc[gene_ID]
        significance_df["padj"]

        for x, padj in zip(significance_df["padj"].index, significance_df["padj"]):
            if padj < 0.05:
                ax.plot(x, bottom - y_range*0.05, color='k', marker=(8, 2, 0), markersize=15, label="Differentially expressed at 5% FDR")
    except KeyError:
        pass

    gene_annotation = annotation_data.loc[gene_ID]

    # title stuff
    ax.set_title(f"{gene_annotation.name}/{gene_annotation['NCBI ID_3']}\n{gene_annotation['Genbank Annotation']}")

@mpl.rc_context({
    'lines.linewidth': 6, 
    'lines.marker':'o', 
    'lines.markersize':18, 
    'legend.fontsize': 'x-large',
    'axes.labelsize': 'xx-large',
    'axes.titlesize':'xx-large',
    'xtick.labelsize':'xx-large',
    'figure.titlesize': 'xx-large',
    'ytick.labelsize':'xx-large'})
def plot_gene_table(gene_df_subset, out_path, rlog_mean_df, rlog_std_df, results_df, annotation_data,
                num_cols = 3,
                night_periods = [(-11, 0), (13, 24)], 
                night_color="#dfdfdf",
                attr_dict={'Control':{'color':'salmon', 'label':'Parental $\it{Prochlorococcus}$'}, 'Pheno':{'color':'lightseagreen', 'label':'Dark-tolerant $\it{Prochlorococcus}$'}}):

    # original from elaina:
    # color_dict={'Control':'#e97e72', 'Pheno':'#52bcc2'}
    if isinstance(gene_df_subset, pd.Series):
        gene_arr = np.array([[gene_df_subset.name]])
    else:
        if len(gene_df_subset) > num_cols:
            gene_arr = list(gene_df_subset.index.values)
            gene_arr += [None]*(num_cols - (len(gene_arr) % num_cols))
            gene_arr = np.array(gene_arr).reshape(-1, num_cols)
        else:
            gene_arr = np.array(gene_df_subset.index.values).reshape(1, len(gene_df_subset))

    y_height = 4
    x_width = 5

    heights = [y_height]*gene_arr.shape[0]
    widths = [x_width]*gene_arr.shape[1]

    fig = plt.figure(figsize=(sum(widths), sum(heights)), constrained_layout=True)
    gs = fig.add_gridspec(ncols=len(widths), nrows=len(heights), height_ratios=heights, width_ratios=widths)

    for i, row in enumerate(gene_arr):
        for j, element in enumerate(row):
            if element != None:
                ax = fig.add_subplot(gs[i,j])
                plot_gene(ax, element, rlog_mean_df, rlog_std_df, results_df, annotation_data, night_periods, night_color, attr_dict)

    # handles, labels = ax.get_legend_handles_labels()
    legend_elements = [mpl.lines.Line2D([0], [0], color=d['color'], label=d['label']) for t, d in attr_dict.items()]
    legend = fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor= (1.01, 0.5))

    xlab = fig.supxlabel("Time (hours)")
    ylab = fig.supylabel("Relative transcript abundance")
    plt.savefig(out_path, bbox_extra_artists=[legend, xlab, ylab], bbox_inches='tight')
    plt.close()

plotting_dir = results_dir / 'all_gene_expression_plots'
plotting_dir.mkdir(exist_ok=True)

# Clock proteins
clock_proteins = {"KaiB":"PMN2A_0914",
                  "KaiC":"PMN2A_0913",
                  "SasA":"PMN2A_0674",
                  "RpaA":"PMN2A_1494",
                  "LdpA":"PMN2A_1131"}

genes = annotation_data[annotation_data["NCBI ID_3"].isin(clock_proteins.values())]
out_path = plotting_dir / 'control_plot.jpg'
plot_gene_table(genes, out_path, rlog_mean_df, rlog_std_df, result_df, annotation_data)

for i, row in annotation_data.iterrows():
    filename = f"{row.name}_{row['NCBI ID_2']}_{row['NCBI ID_3']}_{row['Genbank Annotation']}"
    out_path = plotting_dir / f'{slugify(filename)}.jpeg'
    plot_gene_table(row, out_path, rlog_mean_df, rlog_std_df, result_df, annotation_data)
