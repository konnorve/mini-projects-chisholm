import pandas as pd
from pathlib import Path
import plotly.express as px
import shutil

ref_dir = Path('/nobackup1/chisholmlab/img_proportal/data/derived_data_files')
chosen_references_dir = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/chosen_references_small")
chosen_ref_gff_dir = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/chosen_reference_gffs")

cycog_gene_table = pd.read_table('/nobackup1b/users/chisholmlab/img_proportal/cycogs/CyCOGs-v6.0/cycog-genes.tsv', index_col='gene_iid')
pro_sccg = set(Path('/home/kve/scripts/mini_projects/2022_scope_gradients_HL_adaptation/sccg_ident/pro_sccg.txt').read_text().splitlines())
syn_sccg = set(Path('/home/kve/scripts/mini_projects/2022_scope_gradients_HL_adaptation/sccg_ident/syn_sccg.txt').read_text().splitlines())


organism_table = cycog_gene_table[['taxon_oid', 'organism']].drop_duplicates().set_index('taxon_oid')['organism'].to_dict()

included_genomes = []
cyano_ref_dir = ref_dir / 'fna'
cyano_gff_dir = ref_dir / 'gff'
cyano_refs = list(cyano_ref_dir.iterdir())
for ref in cyano_refs:
    genome_name = ref.stem
    print(genome_name)
    img_id = int(genome_name.split('_IMG_', 1)[1])
    if img_id in organism_table.keys() and organism_table[img_id] in ['Prochlorococcus', 'Synechococcus']:
        cycog_subset = set(cycog_gene_table[cycog_gene_table['taxon_oid']==img_id]['cycog_iid'].to_list())

        if organism_table[img_id] == 'Prochlorococcus':
            pct_sccg = len(cycog_subset.intersection(pro_sccg)) / len(pro_sccg)
        else:
            pct_sccg = len(cycog_subset.intersection(syn_sccg)) / len(syn_sccg)

        print(f"{img_id}\t{pct_sccg}")

        if pct_sccg > 0.1:
            shutil.copy(ref, chosen_references_dir)
            # shutil.copy(cyano_gff_dir / f"{ref.stem}.gff", chosen_ref_gff_dir)
            # included_genomes.append(ref)

# with open('included_pro_syn_genomes.txt', 'w') as f:
#     f.write("\n".join(map(str, included_genomes)))