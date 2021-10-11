import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.colors
from pathlib import Path
from Bio import SeqIO, SeqUtils
import gffpandas.gffpandas as gffpd

# Parsons iMac path:
diel_trained_proj_dir = Path('/nobackup1/kve/2021_trained_diel_expression_vs_genome_location')

trained_diel_df_path = diel_trained_proj_dir / 'proGeneCounts_darkAdaptedControl_vstNormalized.csv'

finalized_df_out_path = diel_trained_proj_dir / 'proGeneCounts_organized_trained_diel_vstNormalized.tsv'

ref_fasta = diel_trained_proj_dir / 'NATL2A_genome_references' / 'onlyNATL2A.fna'
ref_genome_record = SeqIO.read(ref_fasta, 'fasta')
ref_genome_seq = ref_genome_record.seq

gff_file = diel_trained_proj_dir / 'NATL2A_genome_references' / 'onlyNATL2A.gff'
annotation = gffpd.read_gff3(gff_file)
attributes_df = annotation.attributes_to_columns()

figure_out_dir = diel_trained_proj_dir / 'Expression_vs_genome_location'
figure_out_dir.mkdir(exist_ok=True)

trained_diel_df = pd.read_csv(trained_diel_df_path)

# adding important columns to df
trained_diel_df['midpoint'] = (trained_diel_df['start'] + trained_diel_df['end']) / 2

trained_diel_df['rounded_time'] = trained_diel_df['time'].apply(lambda x: int(x))

trained_diel_df['rounded_time_str'] = trained_diel_df['rounded_time'].apply(lambda x: '{:02}'.format(x))

def assign_replicate(sample_num):
    rep_num = (sample_num % 3) + 1
    return rep_num

def assign_sample_num(sample_str):
    return int(sample_str[6:])
    
trained_diel_df['sample_num'] = trained_diel_df['sample'].apply(lambda x: assign_sample_num(x))
trained_diel_df['replicate_num'] = trained_diel_df['sample_num'].apply(lambda x: assign_replicate(x))

trained_diel_df['unqiue_group'] = trained_diel_df['Date'].astype('str') + '_' + trained_diel_df['rounded_time_str'].astype('str') + '_' + trained_diel_df['Treatment']

# make samples dataframe in order to reformat replicates
trained_diel_df_samples = trained_diel_df.copy()
trained_diel_df_samples = trained_diel_df_samples.drop(['gene', 'start', 'end', 'geneFunction', 'sample', 'normTransCounts', 'time', 'midpoint'], axis='columns')
trained_diel_df_samples = trained_diel_df_samples.drop_duplicates()


new_dfs_list = []
for group in trained_diel_df['unqiue_group'].unique():
    df_unique_group = trained_diel_df[trained_diel_df['unqiue_group']==group]
    new_df_slice = df_unique_group.copy()
    new_df_slice = new_df_slice.drop(['normTransCounts', 'replicate_num', 'sample_num', 'sample'], axis='columns')
    new_df_slice = new_df_slice.drop_duplicates()
    new_df_slice = new_df_slice.reset_index()

    unique_replicates = df_unique_group['replicate_num'].unique()

    
    rep_cols = []
    for rep in unique_replicates:
        # print("{} \t {}".format(group, rep))
        rep_col_name = 'rep{}_normTransCounts'.format(rep)
        rep_cols.append(rep_col_name)
        new_df_slice[rep_col_name] = df_unique_group[df_unique_group['replicate_num']==rep].reset_index()['normTransCounts']

    new_df_slice['rep_average_normTransCounts'] = new_df_slice[rep_cols].mean(axis=1)

    new_dfs_list.append(new_df_slice)


organized_trained_diel_df = pd.concat(new_dfs_list)

organized_trained_diel_df = organized_trained_diel_df.set_index(['Treatment', 'unqiue_group', 'gene'])

trained_diel_df.set_index

gene_averages_for_normalization = trained_diel_df.pivot_table(values='normTransCounts', index=['Treatment', 'gene'], aggfunc=np.mean)

gene_averages_for_normalization = gene_averages_for_normalization.rename(columns={'normTransCounts':'normTransCounts_gene_treatment_avg'})

organized_trained_diel_df = organized_trained_diel_df.join(gene_averages_for_normalization)

organized_trained_diel_df = organized_trained_diel_df.reset_index()
organized_trained_diel_df = organized_trained_diel_df.set_index(['Treatment', 'rounded_time', 'gene'])
organized_trained_diel_df = organized_trained_diel_df.sort_values(['Treatment', 'rounded_time', 'gene'])


organized_trained_diel_df['norm_change_normTransCounts'] = organized_trained_diel_df['rep_average_normTransCounts'] / organized_trained_diel_df['normTransCounts_gene_treatment_avg']



def plot_expression_pattern_bar(df_slice, out_path, group_name, width=5000, alpha=0.5, ylims=(0, 2)):
    df_plot = df_slice.sort_values('midpoint')
    
    fig = plt.figure(figsize=(20, 5), constrained_layout=True)
    ax = fig.gca()
    ax.bar(x=df_plot['midpoint'], height=df_plot['norm_change_normTransCounts'], width=width, alpha=alpha)
    ax.set_ylim(bottom=ylims[0], top=ylims[1])

    ax.axhline(y=1, lw=1, c='r')

    ax.set_xlabel('Location on Genome (mb)')
    ax.set_ylabel('NormTransCounts change from treatment average')
    ax.set_title('{} Normalized Gene Counts on Genome Location'.format(group_name))

    plt.savefig(out_path)
    plt.close()

"""
If you could plot as points -- DONE
maybe even connecting adjacent genomic positions and  -- DONE
apply a smoothing function, 

# then (a) you can probably see more clearly any trends that are present and (b) you could fit multiple time points on one graph. 
Chromosome position needs to be converted to polar coordinates/position on a unit circle because it is a circular chromosome.
"""

def smooth_over_count(raw_arr, window_size = 5):
    len_arr = len(raw_arr)
    smoothed_arr = np.zeros(len_arr)
    l_r_len = window_size//2
    for i in range(len_arr):
        indicies_2_smooth = list(range(i - l_r_len, i + l_r_len + 1))
        smoothed_arr[i] = np.mean([raw_arr[j%len_arr] for j in indicies_2_smooth])

    return smoothed_arr

def smooth_over_nucleotides(raw_counts_arr, raw_placement_arr, window_size = 100):

    len_arr = len(raw_counts_arr)
    len_genome = max(raw_placement_arr)
    smoothed_arr = np.zeros(len_arr)
    l_r_len = window_size//2

    assert window_size < len_genome

    for i, i_pos in enumerate(raw_placement_arr):
        
        left_bound = i_pos - l_r_len
        right_bound = i_pos + l_r_len + 1
        counts_2_smooth = []
        if 0 <= left_bound <= len_genome and 0 <= right_bound <= len_genome:
            for pos, count in zip(raw_placement_arr, raw_counts_arr):
                if pos >= left_bound and pos <= right_bound:
                    counts_2_smooth.append(count)

        elif left_bound < 0 and 0 <= right_bound <= len_genome:
            for pos, count in zip(raw_placement_arr, raw_counts_arr):
                if (left_bound % len_genome) <= pos <= len_genome:
                    counts_2_smooth.append(count)
                elif 0 <= pos <= right_bound:
                    counts_2_smooth.append(count)

        elif 0 <= left_bound <= len_genome and right_bound > len_genome:
            for pos, count in zip(raw_placement_arr, raw_counts_arr):
                if left_bound <= pos <= len_genome:
                    counts_2_smooth.append(count)
                elif 0 <= pos <= (right_bound % len_genome):
                    counts_2_smooth.append(count)
        else:
            raise IndexError('right or left values both out of range')

        smoothed_arr[i] = np.mean(counts_2_smooth)

    return smoothed_arr


def plot_expression_pattern(df_treatment, out_path, group_name, window_size = 100, ylims=None):
    
    
    fig = plt.figure(figsize=(20, 5), constrained_layout=True)
    ax = fig.gca()

    unique_groups = df_treatment.index.get_level_values('rounded_time').unique()
    for group in unique_groups:
        df_slice = df_treatment.loc[group]
        df_plot = df_slice.sort_values('midpoint')

        raw_counts = df_plot['norm_change_normTransCounts'].to_numpy()
        raw_positions = df_plot['midpoint'].to_numpy()
        smoothed_counts = smooth_over_nucleotides(raw_counts, raw_positions, window_size=window_size)

        ax.plot(raw_positions, smoothed_counts, label=group)
    
    ax.legend()

    if ylims:
        ax.set_ylim(bottom=ylims[0], top=ylims[1])

    ax.axhline(y=1, lw=1, c='r')

    ax.set_xlabel('Location on Genome (mb)')
    ax.set_ylabel('NormTransCounts change from treatment average')
    ax.set_title('{} Normalized Gene Counts on Genome Location -- window size: {:,}nt'.format(group_name, window_size))

    plt.savefig(out_path)
    plt.close()


def plot_expression_pattern_circle_naive(df_treatment, out_path, group_name, window_size = 100000):
    
    fig = plt.figure(figsize=(20, 20), constrained_layout=True)
    ax = fig.gca()

    max_radius = 0

    unique_groups = df_treatment.index.get_level_values('rounded_time').unique()
    for group in unique_groups:
        df_slice = df_treatment.loc[group]
        df_plot = df_slice.sort_values('midpoint')

        raw_counts = df_plot['norm_change_normTransCounts'].to_numpy()
        raw_positions = df_plot['midpoint'].to_numpy()
        smoothed_counts = smooth_over_nucleotides(raw_counts, raw_positions, window_size=window_size)

        if max(smoothed_counts) > max_radius:
            max_radius = max(smoothed_counts)

        radian_pos = raw_positions * 2 * np.pi / max(raw_positions)

        x_pos = smoothed_counts * np.cos(radian_pos)
        y_pos = smoothed_counts * np.sin(radian_pos)

        ax.plot(x_pos, y_pos, label=group)
    
    ax.legend()

    ax.add_patch(Circle((0,0), radius=1, fill=False))
    
    ax.set_xlim(left=-max_radius, right=max_radius)
    ax.set_ylim(bottom=-max_radius, top=max_radius)

    # ax.axhline(y=1, lw=1, c='r')

    # ax.set_xlabel('Location on Genome (mb)')
    # ax.set_ylabel('NormTransCounts change from treatment average')
    ax.set_title('{} Normalized Gene Counts on Genome Location -- window size: {:,}nt'.format(group_name, window_size))

    plt.savefig(out_path)
    plt.close()


def plot_expression_pattern_circle(df_treatment, out_path, group_name, window_size = 100000, explode = False):
    df_treatment = df_treatment.copy()

    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    plt.axes(projection = 'polar')
    ax = fig.gca()

    max_radius = 0
    genome_length = max(df_treatment['end'])

    unique_groups = df_treatment.index.get_level_values('rounded_time').unique()

    if explode:
        explosion_factor = 0
        for group in unique_groups:
            df_slice = df_treatment.loc[group]
            df_plot = df_slice.sort_values('midpoint')
            raw_counts = df_plot['norm_change_normTransCounts'].to_numpy()
            raw_positions = df_plot['midpoint'].to_numpy()
            smoothed_counts = smooth_over_nucleotides(raw_counts, raw_positions, window_size=window_size)
            if max(smoothed_counts) - min(smoothed_counts) > explosion_factor:
                explosion_factor = max(smoothed_counts) - min(smoothed_counts)

        expansion_factors = [1+(explosion_factor*i) for i in range(len(unique_groups))]
        for group, expansion_factor in zip(unique_groups, expansion_factors):
            df_treatment.loc[group]['norm_change_normTransCounts'] = df_treatment.loc[group]['norm_change_normTransCounts'] * expansion_factor

    
    for group in unique_groups:
        df_slice = df_treatment.loc[group]
        df_plot = df_slice.sort_values('midpoint')

        raw_counts = df_plot['norm_change_normTransCounts'].to_numpy()
        raw_positions = df_plot['midpoint'].to_numpy()
        smoothed_counts = smooth_over_nucleotides(raw_counts, raw_positions, window_size=window_size)

        if max(smoothed_counts) > max_radius:
            max_radius = max(smoothed_counts)

        radian_pos = raw_positions * 2 * np.pi / max(raw_positions)

        # x_pos = smoothed_counts * np.cos(radian_pos)
        # y_pos = smoothed_counts * np.sin(radian_pos)

        ax.plot(radian_pos, smoothed_counts, label=group)
    
    ax.legend()

    xticks = []
    x = 0
    while x < genome_length:
        xticks.append(x)
        x+=250000
    xticks_rad = [x*2*np.pi/genome_length for x in xticks]
    ax.set_xticks(xticks_rad)
    ax.set_xticklabels(["{:,}".format(x) for x in xticks])

    if explode:
        ax.set_yticks(expansion_factors)
        ax.set_yticklabels(unique_groups)

    # ax.set_xlabel('Location on Genome (mb)')
    # ax.set_ylabel('NormTransCounts change from treatment average')
    ax.set_title('{} Normalized Gene Counts on Genome Location -- window size: {:,}nt'.format(group_name, window_size))

    plt.savefig(out_path)
    plt.close()

def plot_expression_pattern_circle_time_exp(df_treatment, out_path, group_name, min_max, window_size = 100000):
    
    df_treatment_copy = df_treatment.copy()
    df_treatment_copy = df_treatment_copy.reset_index()
    dfs_to_plot = []
    if min_max == 'max':
        dfs_to_plot.append(df_treatment_copy.sort_values('norm_change_normTransCounts', ascending=False).drop_duplicates(['gene']))
    elif min_max == 'min':
        dfs_to_plot.append(df_treatment_copy.sort_values('norm_change_normTransCounts', ascending=True).drop_duplicates(['gene']))
    elif min_max == 'both':
        dfs_to_plot.append(df_treatment_copy.sort_values('norm_change_normTransCounts', ascending=False).drop_duplicates(['gene']))
        dfs_to_plot.append(df_treatment_copy.sort_values('norm_change_normTransCounts', ascending=True).drop_duplicates(['gene']))

    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    plt.axes(projection = 'polar')
    ax = fig.gca()


    for df_treatment_copy in dfs_to_plot:
        max_radius = 0
        genome_length = max(df_treatment_copy['end'])

        df_plot = df_treatment_copy.sort_values('midpoint')

        raw_counts = df_plot['norm_change_normTransCounts'].to_numpy()
        raw_positions = df_plot['midpoint'].to_numpy()
        times = df_plot['rounded_time'].to_numpy()
        smoothed_counts = smooth_over_nucleotides(raw_counts, raw_positions, window_size=window_size)
        
        unique_times = np.unique(times)

        cmap = plt.cm.viridis
        norm = matplotlib.colors.Normalize(vmin=min(unique_times), vmax=max(unique_times))

        if max(smoothed_counts) > max_radius:
            max_radius = max(smoothed_counts)

        radian_pos = raw_positions * 2 * np.pi / max(raw_positions)

        for t in unique_times:
            ix = np.where(times == t)
            ax.scatter(radian_pos[ix], smoothed_counts[ix], color=cmap(norm(t)), label=t, alpha=0.5)
    
    ax.legend()

    xticks = []
    x = 0
    while x < genome_length:
        xticks.append(x)
        x+=250000
    xticks_rad = [x*2*np.pi/genome_length for x in xticks]
    ax.set_xticks(xticks_rad)
    ax.set_xticklabels(["{:,}".format(x) for x in xticks])

    # ax.set_xlabel('Location on Genome (mb)')
    # ax.set_ylabel('NormTransCounts change from treatment average')
    ax.set_title('{} Normalized Gene Counts on Genome Location -- window size: {:,}nt'.format(group_name, window_size))

    plt.savefig(out_path)
    plt.close()

def plot_expression_pattern_circle_min_max(df_treatment, out_path, group_name, window_size = 100000):
    
    df_treatment_copy = df_treatment.copy()
    df_treatment_copy = df_treatment_copy.reset_index()

    df_treatment_copy_max = df_treatment_copy.sort_values('norm_change_normTransCounts', ascending=False).drop_duplicates(['gene'])
    df_treatment_copy_min = df_treatment_copy.sort_values('norm_change_normTransCounts', ascending=True).drop_duplicates(['gene'])

    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    plt.axes(projection = 'polar')
    ax = fig.gca()

    for df_treatment_copy in [df_treatment_copy_max, df_treatment_copy_min]:

        max_radius = 0
        genome_length = max(df_treatment_copy['end'])

        df_plot = df_treatment_copy.sort_values('midpoint')

        raw_counts = df_plot['norm_change_normTransCounts'].to_numpy()
        raw_positions = df_plot['midpoint'].to_numpy()
        times = df_plot['rounded_time'].to_numpy()
        smoothed_counts = smooth_over_nucleotides(raw_counts, raw_positions, window_size=window_size)
        
        unique_times = np.unique(times)

        cmap = plt.cm.viridis
        norm = matplotlib.colors.Normalize(vmin=min(unique_times), vmax=max(unique_times))

        if max(smoothed_counts) > max_radius:
            max_radius = max(smoothed_counts)

        radian_pos = raw_positions * 2 * np.pi / max(raw_positions)

        for t in unique_times:
            ix = np.where(times == t)
            ax.scatter(radian_pos[ix], smoothed_counts[ix], color=cmap(norm(t)), label=t, alpha=0.5)
    
        ax.legend()

    xticks = []
    x = 0
    while x < genome_length:
        xticks.append(x)
        x+=250000
    xticks_rad = [x*2*np.pi/genome_length for x in xticks]
    ax.set_xticks(xticks_rad)
    ax.set_xticklabels(["{:,}".format(x) for x in xticks])

    # ax.set_xlabel('Location on Genome (mb)')
    # ax.set_ylabel('NormTransCounts change from treatment average')
    ax.set_title('{} Normalized Gene Counts on Genome Location -- window size: {:,}nt'.format(group_name, window_size))

    plt.savefig(out_path)
    plt.close()


def get_GC_arr(ref_seq, window_size):

    genome_len = len(ref_seq)
    gc_smooth_arr = np.zeros(genome_len)
    l_r_len = window_size // 2

    def convert2gc(b):
        if b in 'ATat':
            return 0
        return 1

    gc_arr = np.array([convert2gc(b) for b in str(ref_seq)])

    assert window_size < genome_len

    for i in range(genome_len):
        left_bound = i - l_r_len
        right_bound = i + l_r_len + 1

        if 0 <= left_bound <= genome_len and 0 <= right_bound < genome_len:
            s_val = np.sum(gc_arr[left_bound:right_bound])
        else:
             s_val = np.sum(gc_arr[left_bound:]) + np.sum(gc_arr[:right_bound%genome_len])

        gc_smooth_arr[i] = s_val / window_size

    return gc_smooth_arr

def plot_GC_content(ref_seq, out_path, window_size, gc_arr=None):
    if gc_arr is None:
        gc_arr = get_GC_arr(ref_seq, window_size)
    global_gc_pct = SeqUtils.GC(ref_seq) / 100
    genome_len = len(ref_seq)
    
    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    plt.axes(projection = 'polar')
    ax = fig.gca()

    radian_pos = np.linspace(0, 2*np.pi, genome_len)

    ax.plot(radian_pos, gc_arr, color='b')

    ax.axvline(radian_pos[gc_arr.argmax()], color='g')
    ax.axvline(radian_pos[gc_arr.argmin()], color='r')
    ax.axhline(global_gc_pct, color='b')

    xticks = []
    x = 0
    while x < genome_len:
        xticks.append(x)
        x+=250000
    xticks_rad = [x*2*np.pi/genome_len for x in xticks]
    ax.set_xticks(xticks_rad)
    ax.set_xticklabels(["{:,}".format(x) for x in xticks])

    ax.set_title('GC Content Averaged over window size: {:,}nt'.format(window_size))

    plt.savefig(out_path)
    plt.close()


def plot_GC_content_histogram(ref_seq, out_path, window_size, gc_arr=None):
    if gc_arr is None:
        gc_arr = get_GC_arr(ref_seq, window_size)
    global_gc_pct = SeqUtils.GC(ref_seq) / 100
    genome_len = len(ref_seq)
    
    fig = plt.figure(figsize=(5, 5), constrained_layout=True)
    ax = fig.gca()

    ax.hist(gc_arr, color='b')
    ax.axvline(global_gc_pct, color='black')

    ax.set_title('GC Content Histogram with Window Size: {:,}nt'.format(window_size))

    plt.savefig(out_path)
    plt.close()


def get_GC_skew_arr(ref_seq, window_size):

    # https://doi.org/10.1093/nar/26.10.2286
    # skew calculated as (G-C)/(G+C)

    genome_len = len(ref_seq)
    gc_smooth_arr = np.zeros(genome_len)
    l_r_len = window_size // 2

    def convert2gcskew(b):
        if b in 'ATat':
            return 0
        elif b in 'Gg':
            return 1
        else:
            return -1
    
    gc_skew_arr = np.array([convert2gcskew(b) for b in str(ref_seq)])
    gc_arr = np.absolute(gc_skew_arr)

    assert window_size < genome_len

    for i in range(genome_len):
        left_bound = i - l_r_len
        right_bound = i + l_r_len + 1

        if 0 <= left_bound <= genome_len and 0 <= right_bound < genome_len:
            gc_smooth_arr[i] = np.sum(gc_skew_arr[left_bound:right_bound]) / np.sum(gc_arr[left_bound:right_bound])
        else:
             gc_smooth_arr[i] = (np.sum(gc_skew_arr[left_bound:]) + np.sum(gc_skew_arr[:right_bound%genome_len])) / (np.sum(gc_arr[left_bound:]) + np.sum(gc_arr[:right_bound%genome_len]))

        if i % 100000 == 0: print(i)

    return gc_smooth_arr


def plot_skew_content(ref_seq, out_path, window_size, gc_skew_arr=None, gc_cum_skew_arr=None, polar=True):

    if gc_skew_arr is None:
        gc_skew_arr = get_GC_skew_arr(ref_seq, window_size)
    if gc_cum_skew_arr is None:
        gc_cum_skew_arr = np.cumsum(gc_skew_arr)

    genome_len = len(ref_seq)
    
    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    if polar:
        plt.axes(projection = 'polar')
        radian_pos = np.linspace(0, 2*np.pi, genome_len)
        ax = fig.gca()
    else:
        ax = fig.gca()
        radian_pos = range(genome_len)
        ax_secondary = ax.twinx()
        ax_secondary.plot(radian_pos, gc_skew_arr, color='blue', label='GC skew', alpha=0.4)
        ax_secondary.axhline(y=0, c='black')

    ax.plot(radian_pos, gc_cum_skew_arr, color='orange', label='Cumulative GC skew')

    ax.axvline(radian_pos[gc_cum_skew_arr.argmax()], color='green')
    ax.axvline(radian_pos[gc_cum_skew_arr.argmin()], color='red')

    ax.legend()

    if polar:
        xticks = []
        x = 0
        while x < genome_len:
            xticks.append(x)
            x+=250000
        xticks_rad = [x*2*np.pi/genome_len for x in xticks]
        ax.set_xticks(xticks_rad)
        ax.set_xticklabels(["{:,}".format(x) for x in xticks])

    ax.set_title('GC Skew and Cumulative Sum over window size: {:,}nt'.format(window_size))

    plt.savefig(out_path)
    plt.close()



window_sizes = [1000, 5000, 10000, 50000, 100000, 500000]

for ws in window_sizes:
    print(ws)

    gc_arr = get_GC_arr(ref_genome_seq, ws)

    (figure_out_dir / 'GC_content').mkdir(exist_ok=True)
    out_path = figure_out_dir / 'GC_content' / 'GC_Content_window_{:06}'.format(ws)
    plot_GC_content(ref_genome_seq, out_path, ws, gc_arr)

    (figure_out_dir / 'GC_content_histograms').mkdir(exist_ok=True)
    out_path = figure_out_dir / 'GC_content_histograms' / 'GC_Content_histogram_window_{:06}'.format(ws)
    plot_GC_content_histogram(ref_genome_seq, out_path, ws, gc_arr=None)

    gc_skew_arr = get_GC_skew_arr(ref_genome_seq, ws)
    gc_cum_skew_arr = np.cumsum(gc_skew_arr)

    (figure_out_dir / 'GC_skew').mkdir(exist_ok=True)
    out_path = figure_out_dir / 'GC_skew' / 'GC_Skew_polar_window_{:06}'.format(ws)
    plot_skew_content(ref_genome_seq, out_path, ws, gc_skew_arr, gc_cum_skew_arr, polar=True)

    (figure_out_dir / 'GC_skew').mkdir(exist_ok=True)
    out_path = figure_out_dir / 'GC_skew' / 'GC_Skew_linear_window_{:06}'.format(ws)
    plot_skew_content(ref_genome_seq, out_path, ws, gc_skew_arr, gc_cum_skew_arr, polar=False)



organized_trained_diel_df.to_csv(finalized_df_out_path, sep='\t')

unique_treatments = organized_trained_diel_df.index.get_level_values('Treatment').unique()

for treatment in unique_treatments:

    trained_diel_df_treatment = organized_trained_diel_df.loc[treatment]
    

    for ws in window_sizes:
        print(treatment, '\t', ws)

        (figure_out_dir / 'expression_pattern_line').mkdir(exist_ok=True)
        out_path = figure_out_dir / 'expression_pattern_line' / '{}_{:06}_expression_pattern_line'.format(treatment, ws)
        plot_expression_pattern(trained_diel_df_treatment, out_path, treatment, window_size=ws)

        (figure_out_dir / 'expression_pattern_circular').mkdir(exist_ok=True)
        out_path = figure_out_dir / 'expression_pattern_circular' / '{}_{:06}_expression_pattern_circular'.format(treatment, ws)
        plot_expression_pattern_circle(trained_diel_df_treatment, out_path, treatment, window_size = ws)

        (figure_out_dir / 'expression_pattern_circular_exploded').mkdir(exist_ok=True)
        out_path = figure_out_dir / 'expression_pattern_circular_exploded' / '{}_{:06}_expression_pattern_circular_exploded'.format(treatment, ws)
        plot_expression_pattern_circle(trained_diel_df_treatment, out_path, treatment, window_size = ws, explode=True)

        (figure_out_dir / 'time_max_exp_circular').mkdir(exist_ok=True)
        out_path = figure_out_dir / 'time_max_exp_circular' / '{}_{:06}_time_max_exp_circular'.format(treatment, ws)
        plot_expression_pattern_circle_time_exp(trained_diel_df_treatment, out_path, treatment, min_max='max', window_size = ws)

        (figure_out_dir / 'time_min_exp_circular').mkdir(exist_ok=True)
        out_path = figure_out_dir / 'time_min_exp_circular' / '{}_{:06}_time_min_exp_circular'.format(treatment, ws)
        plot_expression_pattern_circle_time_exp(trained_diel_df_treatment, out_path, treatment, min_max='min', window_size = ws)

        (figure_out_dir / 'time_min_max_exp_circular').mkdir(exist_ok=True)
        out_path = figure_out_dir / 'time_min_max_exp_circular' / '{}_{:06}_time_min_max_exp_circular'.format(treatment, ws)
        plot_expression_pattern_circle_time_exp(trained_diel_df_treatment, out_path, treatment, min_max='both', window_size = ws)