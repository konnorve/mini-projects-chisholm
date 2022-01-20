from pathlib import Path
import pandas as pd
import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

results_dir = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/results_eggNOG_ncbi_gff_10_mapq")

experiment_4_results_df = pd.read_csv(results_dir / "experiments" / "experiment_4" / "DGE_tables" / "experiment_4_MIT9301_DGE_all.tsv", index_col=0, sep="\t")
experiment_11_results_df = pd.read_csv(results_dir / "experiments" / "experiment_11" / "DGE_tables" / "experiment_11_MIT9301_DGE_all.tsv", index_col=0,  sep="\t")

meta_columns = ['product', 'seq_id', 'source', 'type', 'start', 'end',
    'score', 'strand', 'phase', 'attributes', 'BRITE', 'BiGG_Reaction',
    'CAZy', 'COG_category', 'Dbxref', 'Description', 'EC', 'GOs',
    'Is_circular', 'KEGG_Module', 'KEGG_Pathway', 'KEGG_Reaction',
    'KEGG_TC', 'KEGG_ko', 'KEGG_rclass', 'Name', 'Note', 'PFAMs', 'Parent',
    'Preferred_name', 'anticodon', 'bound_moiety', 'eggNOG_OGs',
    'end_range', 'evalue', 'exception', 'gbkey', 'gene', 'gene_biotype',
    'genome', 'inference', 'locus_tag', 'max_annot_lvl', 'mol_type',
    'old_locus_tag', 'partial', 'protein_id', 'pseudo', 'regulatory_class',
    'score_annot', 'seed_ortholog', 'strain', 'transl_table', 'organism',
    'gene_synonym', 'start_range']

concat_df_inner = experiment_11_results_df.merge(experiment_4_results_df, how='inner', on=meta_columns, suffixes=('_exp11', '_exp4'))
concat_df_outer = experiment_11_results_df.merge(experiment_4_results_df, how='inner', on=meta_columns, suffixes=('_exp11', '_exp4'))
assert concat_df_inner.shape == concat_df_outer.shape

concat_df = concat_df_inner

col_print_list = '\n'.join([str(col) for col in concat_df.columns.to_list()])
log.debug(f"concat_df columns:\n{col_print_list}")
log.debug(f"concat_df:\n{concat_df}")

concat_df.to_csv("concat_df.tsv", sep='\t')