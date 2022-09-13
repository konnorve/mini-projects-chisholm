import pandas as pd
from pathlib import Path
import math
from itertools import combinations
import matplotlib.pyplot as plt
from venn import venn

df_dict = {
    "bcftools_all" : pd.read_table("/nfs/chisholmlab001/kve/2022_WH7803_variant_calling/og_vs_ax_x_results/bcftools_all/all_genomes_filtered_variants.bcftools_all.tsv"),
    "bcftools_std" : pd.read_table("/nfs/chisholmlab001/kve/2022_WH7803_variant_calling/og_vs_ax_x_results/bcftools_standard/all_genomes_filtered_variants.bcftools_standard.tsv"),
    "freebayes" : pd.read_table("/nfs/chisholmlab001/kve/2022_WH7803_variant_calling/og_vs_ax_x_results/freebayes/all_genomes_filtered_variants.freebayes.tsv"),
}

STD_VCF_COLUMNS = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']

for k in df_dict.keys():
    df_dict[k]['ident'] = df_dict[k].apply(lambda r: "|".join([str(r[col]) for col in ['POS', 'REF', 'ALT']]), axis=1)
    df_dict[k] = df_dict[k].set_index('ident')

def venn_set(include_sets, exclude_sets, set_dict):
    included = set.intersection(*[set_dict[i] for i in include_sets])
    excluded = set.union(*[set_dict[e] for e in exclude_sets])
    return included.difference(excluded)

def venn_name(include_sets, exclude_sets, sep=', '):
    return sep.join([f"+{i}" for i in include_sets]+[f"-{e}" for e in exclude_sets])

def sets_to_venn(set_dict):
    venn_dict = dict()
    set_keys = set(set_dict.keys())
    for i in range(1,len(set_keys)+1):
        for group in combinations(set_keys, i):
            include = set(group)
            exclude = set_keys-set(group)
            if len(exclude) > 0:
                venn_dict[venn_name(include, exclude)] = venn_set(include, exclude, set_dict)
            else:
                name = ", ".join([f"+{i}" for i in include])
                venn_dict[name] = set.intersection(*[set_dict[i] for i in include])
    return venn_dict

samples = set([i for df in df_dict.values() for i in df.genome_name.to_list()])
for sample in samples:
    print(sample)
    set_dict = {k : set(df[df.genome_name == sample].index.to_list()) for k,df in df_dict.items()}
    venn_dict = sets_to_venn(set_dict)

    for k,s in set_dict.items():
        print(k, len(s), sep='\t')
    for k,s in venn_dict.items():
        print(k, len(s), sep='\t')
        for n,df in df_dict.items():
            if f"+{n}" in k:
                df[df.genome_name==sample][STD_VCF_COLUMNS].rename(columns={'CHROM':'#CHROM'}).loc[s].to_csv(f"sets/{sample}_{k.replace(', ','')}.vcf", sep='\t', index=False)
                break
    
    fig = plt.figure()
    venn(set_dict, ax=fig.gca())
    plt.savefig(f"{sample}_venn.png")

