import pandas as pd
from pathlib import Path
import numpy as np



scratch_dir = Path('/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2')

flagstat_counts = pd.read_table("read_count_results.tsv", header=[0,1], index_col=0)
flagstat_counts = flagstat_counts[['read_counts','big_index']]
gene_2_cycog = pd.read_table("/nobackup1b/users/chisholmlab/img_proportal/cycogs/CyCOGs-v6.0/cycog-genes.tsv", index_col='gene_iid')['cycog_iid'].to_dict()

for sample in (scratch_dir / '2_htseq_output').iterdir():
    if sample.is_dir():
        print(sample.name)
        cycog_count_series = []
        metadata_series = []
        for genome in sample.iterdir():
            if genome.suffix == '.tsv':
                genome_counts_df = pd.read_table(genome, header=0, names=['gene_id', 'paired_counts', 'unpaired_counts'])
                genome_metadata_df = genome_counts_df[genome_counts_df['gene_id'].str.contains("__")]
                genome_counts_df = genome_counts_df[~genome_counts_df['gene_id'].str.contains("__")]
                genome_counts_df['gene_id'] = genome_counts_df['gene_id'].map(gene_2_cycog.get)
                genome_counts_df = genome_counts_df.groupby("gene_id").sum()

                genome_metadata_df['gene_id'] = genome_metadata_df['gene_id'].apply(lambda x: x[2:])
                genome_metadata_series = genome_metadata_df.set_index('gene_id').sum(axis=1)
                genome_metadata_series.name = genome.stem

                genome_counts_series = genome_counts_df.sum(axis=1)
                genome_counts_series.name = genome.stem
                cycog_count_series.append(genome_counts_series)
                metadata_series.append(genome_metadata_series)
        df = pd.concat(cycog_count_series, axis='columns')
        df = df.replace(np.nan, 0)
        print(df)
        print(df.sum().sum())

        df = pd.concat(metadata_series, axis='columns')
        print(df)
        break
