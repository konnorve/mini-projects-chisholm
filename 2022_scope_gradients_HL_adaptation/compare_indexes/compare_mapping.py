import pandas as pd
from pathlib import Path
import plotly.express as px

read_counts_path = Path("read_counts.tsv")
big_index_flagstat_dir = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/3_flagstat")
small_index_flagstat_dir = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/3_flagstat_small_refset")

read_counts = pd.read_table(read_counts_path, names=['sample', 'type', 'count'])
read_counts = read_counts.pivot(index="sample", columns="type", values="count")

def get_df(d):
    df_arr = []
    for path in d.iterdir():
        if path.suffix == '.tsv':
            count = path.read_text().split('\t', 1)[0]
            stem = path.stem.split("_")
            if len(stem)==2:
                df_arr.append([f"{stem[0]}_{stem[1]}", "all", count])
            elif len(stem)==3:
                if stem[2] == 'paired':
                    df_arr.append([f"{stem[0]}_{stem[1]}", "paired", count])
                elif stem[2] == 'unpaired':
                    df_arr.append([f"{stem[0]}_{stem[1]}", "unpaired", count])
                else:
                    raise ValueError
            else:
                raise ValueError
    df = pd.DataFrame(df_arr, columns=['sample', 'type', 'count'])
    df = df.pivot(index="sample", columns="type", values="count")
    return df


d = {
    'read_counts' : read_counts,
    'big_index' : get_df(big_index_flagstat_dir), 
    'small_index' : get_df(small_index_flagstat_dir),
}
df = pd.concat(d.values(), axis=1, keys=d.keys())
df = df.astype('float')

df['read_counts', 'all'] = df.loc[:,'read_counts'].sum(axis=1)
df['analysis','pct_reads_mapped_big_index'] = df.loc[:,('big_index', 'all')] / df.loc[:,('read_counts', 'all')]
df['analysis','pct_reads_mapped_small_index'] = df.loc[:,('small_index', 'all')] / df.loc[:,('read_counts', 'all')]
df['analysis','pct_mapped_reads_small_of_big_index'] = df.loc[:,('small_index', 'all')] / df.loc[:,('big_index', 'all')]

px.histogram(df['analysis'], x="pct_reads_mapped_big_index").write_image("pct_reads_mapped_big_index.png")
px.histogram(df['analysis'], x="pct_reads_mapped_small_index").write_image("pct_reads_mapped_small_index.png")
px.histogram(df['analysis'], x="pct_mapped_reads_small_of_big_index").write_image("mapped_reads_ratio_small_vs_big_index.png")

df.to_csv("read_count_results.tsv", sep="\t")

df = df.dropna()

print(df['analysis','pct_reads_mapped_small_index'].sum())
print(df['analysis','pct_reads_mapped_big_index'].sum())