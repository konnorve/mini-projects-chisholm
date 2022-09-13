import pandas as pd
import numpy as np
from pathlib import Path
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def agg_vcf_df(files):
    dfs2concat = []
    for sample_path in files:
        sample_name = sample_path.parent.parent.name
        logging.info(sample_name)
        sample_df = pd.read_csv(sample_path, sep='\t', comment='#', names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
        logging.debug(sample_df)
        col_order = list(sample_df.columns)
        sample_df['genome_name'] = sample_name
        col_order.insert(0, 'genome_name')
        sample_df = sample_df[col_order]
        logging.info(len(sample_df))
        logging.debug(sample_df)
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


vcf_dir = Path('/nfs/chisholmlab001/sbiller/archived_projects/allison_dark_phenotype')
files = sorted(list(Path(vcf_dir).glob("*/output/output.vcf")))
print("\n".join(map(str, files)))
vcf_df = agg_vcf_df(files)
print(vcf_df['genome_name'].unique())
vcf_df = vcf_df[vcf_df['ALT']!='.']
vcf_df = convert_dict_column_to_columns(vcf_df, 'INFO', ';', '=')


vcf_df = vcf_df[vcf_df['genome_name'].isin(
    ['7932_poly_mutref_output_v352',
    '7933_poly_mutref_output_v352',
    '7934_poly_mutref_output_v352',
    '7935_poly_mutref_output_v352',
    '7936_poly_mutref_output_v352',
    '7937_poly_mutref_output_v352']
)]

vcf_df = vcf_df.replace({
    '7932_poly_mutref_output_v352':"N2_A1A_before_dark_A",
    '7933_poly_mutref_output_v352':"N2_A1A_before_dark_B",
    '7934_poly_mutref_output_v352':"N2_A1A_afterT1_A",
    '7935_poly_mutref_output_v352':"N2_A1A_afterT1_B",
    '7936_poly_mutref_output_v352':"N2_A1A_afterT6_A",
    '7937_poly_mutref_output_v352':"N2_A1A_afterT6_B"
})

vcf_df.to_csv('/nfs/chisholmlab001/kve/2022_SNPs_Dark_Adapted_Genomes/analysis/steve_variant_calls.breseq.tsv', sep='\t', index=False)