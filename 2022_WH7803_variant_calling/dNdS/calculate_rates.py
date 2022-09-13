import pandas as pd
import numpy as np
from pathlib import Path
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def agg_vcf_df(files):
    dfs2concat = []
    for sample_path in files:
        sample_name = "_".join(Path(sample_path).name.split('.')[0:2])
        logging.info(sample_name)
        sample_df = pd.read_csv(sample_path, sep='\t', comment='#', header=None, usecols=range(9))
        sample_df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        logging.debug(sample_df)
        col_order = list(sample_df.columns)
        sample_df['genome_name'] = sample_name
        col_order.insert(0, 'genome_name')
        sample_df = sample_df[col_order]
        logging.info(sample_df)
        dfs2concat.append(sample_df)
    return pd.concat(dfs2concat, axis=0, ignore_index=True)

def convert_dict_column_to_columns(df, col, element_deliminator, key_value_deliminator):
    
    out_col_names = get_unique_INFO_elements(df)

    col_list = df[col].to_numpy()
    col_arr = np.empty((len(col_list), len(out_col_names)), dtype=object)
    for i in range(len(col_list)):
        row_dict = dict(ele.split("=") for ele in col_list[i].split(element_deliminator) if len(ele.split(key_value_deliminator))==2)
        for j, col in enumerate(out_col_names):
            col_arr[i][j] = row_dict.get(col)

    df[out_col_names] = col_arr

    return df

def get_unique_INFO_elements(df):
    info_list = df['INFO'].to_numpy()
    
    info_columns = []
    for i in range(info_list.size):
        s = info_list[i]
        for ele in s.split(';'):
            col = ele[:ele.find('=')]
            if col not in info_columns:
                info_columns.append(col)

    return info_columns

def dNdS(annotated_vcf_dir, dNdS_table_outpath):
    files = list(Path(annotated_vcf_dir).glob("*.vcf"))
    vcf_df = agg_vcf_df(files)
    vcf_df = vcf_df[vcf_df['ALT']!='.']
    vcf_df = convert_dict_column_to_columns(vcf_df, 'INFO', ';', '=')
    vcf_df = vcf_df.groupby(['genome_name','IsSynonymous']).size().unstack(fill_value=0)
    vcf_df = vcf_df.rename(columns={'0':'nonsynonymous', '1':'synonymous','9':'N/A or Unknown'})
    vcf_df['dN/dS'] = vcf_df['nonsynonymous'] / vcf_df['synonymous']
    vcf_df.to_csv(dNdS_table_outpath, sep='\t')

dNdS("/nfs/chisholmlab001/kve/2022_WH7803_variant_calling/og_vs_ax_x_results/annotated_vcfs", 'og_vs_ax_x_dNdS.tsv')
dNdS("/nfs/chisholmlab001/kve/2022_WH7803_variant_calling/ax_vs_x_results/annotated_vcfs", 'ax_vs_x_dNdS.tsv')
