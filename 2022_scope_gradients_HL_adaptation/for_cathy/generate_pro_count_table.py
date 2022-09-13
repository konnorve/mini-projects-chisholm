import pandas as pd
from pathlib import Path

img_genome_df = pd.read_table("/nobackup1b/users/chisholmlab/img_proportal/data/lists/img_genome_lookup_table.txt")
cycog_genes = pd.read_table("/nobackup1b/users/chisholmlab/img_proportal/cycogs/CyCOGs-v6.0/cycog-genes.tsv", index_col='gene_iid')

pro_sccg_cycog_list = Path("/home/kve/scripts/mini_projects/2022_scope_gradients_HL_adaptation/sccg_ident/pro_sccg.txt").read_text().split("\n")

htseq_dir = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/2_htseq_concat_output_0_mapq")

htseq_metadata_list = []
htseq_count_list = []
for i, f in enumerate(Path(htseq_dir).iterdir()):
    print(i, f)
    if f.suffix == '.tsv':
        htseq_df = pd.read_table(f, header=0, names=['gene_id', 'paired_counts', 'unpaired_counts'])
        htseq_df.loc[:,'sample'] = f.stem

        htseq_counts_df = htseq_df[~htseq_df['gene_id'].str.contains("__")].copy()
        htseq_counts_df.loc[:,'counts'] = htseq_counts_df['paired_counts'] + htseq_counts_df['unpaired_counts']
        htseq_counts_df = htseq_counts_df.join(cycog_genes, on='gene_id')
        htseq_counts_df = htseq_counts_df.set_index(['sample','cycog_iid','gene_id'])

        htseq_metadata_df = htseq_df[htseq_df['gene_id'].str.contains("__")].copy()
        htseq_metadata_df.loc[:,'gene_id'] = htseq_metadata_df['gene_id'].str.replace("__", "")
        htseq_metadata_df = htseq_metadata_df.set_index(['sample','gene_id'])
        htseq_metadata_df.loc[(f.stem, 'read_counts'), :] = htseq_counts_df.sum(axis=0)

        htseq_metadata_list.append(htseq_metadata_df)
        htseq_count_list.append(htseq_counts_df)

htseq_metadata_df = pd.concat(htseq_metadata_list)
htseq_counts_df = pd.concat(htseq_count_list)

htseq_counts_df = htseq_counts_df[htseq_counts_df['organism']=='Prochlorococcus']

htseq_counts_df = htseq_counts_df.groupby(level=['sample','cycog_iid'])['counts'].sum().unstack(level=0)
sccg_avg_series = htseq_counts_df.loc[pro_sccg_cycog_list].mean(axis=0)
htseq_counts_df = htseq_counts_df.divide(sccg_avg_series, axis='columns')

htseq_counts_df.to_csv('prochlorococcus_cycog_count_table.csv')