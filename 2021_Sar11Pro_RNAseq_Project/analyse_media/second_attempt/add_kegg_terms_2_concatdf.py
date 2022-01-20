from pathlib import Path
import pandas as pd
import numpy as np
import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

concat_df = pd.read_csv("concat_df.tsv", sep='\t')
ref_table = pd.read_csv("kegg_ids.tsv", sep='\t')

['ref_col', 'kegg_df_ref', 'kegg_id', 'kegg_term']

concat_df_kegg_cols = ref_table.ref_col.unique()

concat_df = concat_df.fillna(value={col:"" for col in concat_df_kegg_cols})
concat_df = concat_df.replace(to_replace="-", value={col:"" for col in concat_df_kegg_cols})

for col in ref_table.ref_col.unique():
    ref_table_subset = ref_table[ref_table['ref_col']==col]
    ref_dict = {row.loc['kegg_id']:row.loc['kegg_term'] for i, row in ref_table_subset.iterrows()}
    # log.debug(f"ref_dict:\n{ref_dict}")
    col_ids = concat_df[col].to_list()
    col_terms = []
    for row in col_ids:
        if row == "":
            row_label = ""
        else:
            row_label = ", ".join([str(ref_dict[i]) for i in str(row).split(',')])
        col_terms.append(row_label)

    concat_df[f"{col}_terms"] = col_terms

concat_df.to_csv("concat_df_w_terms.tsv", sep='\t')
