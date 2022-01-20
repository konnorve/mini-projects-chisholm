from pathlib import Path
import pandas as pd
import numpy as np
import re
import logging as log
import urllib.error
import pickle

import Bio.KEGG.REST as kegg_rest

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

concat_df = pd.read_csv("concat_df.tsv", sep='\t')

kegg_pathway_columns = {'BRITE':'brite', 'KEGG_Module':'module', 'KEGG_Pathway':'pathway', 'KEGG_Reaction':'reaction', 'KEGG_rclass':'rclass'}

def ident_kegg_term(db, id):
    try: 
        return kegg_rest.kegg_find(db, id).readline().split('\n')[0].split('\t')[1]
    except IndexError as e:
        log.error(f'problematic term:\n{kegg_rest.kegg_find(db, id).readline()}')
        return ""
    except urllib.error.URLError as e:
        return e.read().decode("utf8", 'ignore')

dfs = []
for col, db in kegg_pathway_columns.items():
    log.info(f"{col}\t{db}")
    unique_terms = set()
    for row in concat_df[col].to_list():
        unique_terms |= set(str(row).split(","))
    unique_terms -= {'-', 'nan'}
    term_table = []
    for i in unique_terms:
        term = ident_kegg_term(db, i)
        log.info(f"{col}\t{i}\t{term}")
        term_table.append([col, db, i, term])
    dfs.append(pd.DataFrame(term_table, columns=['ref_col', 'kegg_df_ref', 'kegg_id', 'kegg_term']))

ref_table = pd.concat(dfs, axis=0)
ref_table.to_csv("kegg_ids.tsv", sep='\t')


