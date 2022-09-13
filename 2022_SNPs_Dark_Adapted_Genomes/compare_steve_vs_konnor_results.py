import pandas as pd

steve_df = pd.read_table("/nfs/chisholmlab001/kve/2022_SNPs_Dark_Adapted_Genomes/analysis/steve_variant_calls.breseq.tsv", sep='\t')
# freebayes
konnor_df = pd.read_table("/nfs/chisholmlab001/kve/2022_SNPs_Dark_Adapted_Genomes/illumina_results/freebayes/all_genomes_filtered_variants.freebayes.tsv", sep='\t')

# # bcftools standard
# konnor_df = pd.read_table("/nfs/chisholmlab001/kve/2022_SNPs_Dark_Adapted_Genomes/analysis/all_genomes_pacbio_illumina.bcftools_standard.tsv", sep='\t')

# # bcftools all
# konnor_df = pd.read_table("/nfs/chisholmlab001/kve/2022_SNPs_Dark_Adapted_Genomes/illumina_results/bcftools_all/all_genomes_filtered_variants.bcftools_all.tsv", sep='\t')

steve_df = steve_df[steve_df['CHROM']=='NC_007335']

illumina_samples = steve_df['genome_name'].unique()

konnor_df = konnor_df[konnor_df['genome_name'].isin(illumina_samples)]

steve_df['ident'] = steve_df.apply(lambda r: f"{r['POS']}|{r['REF']}|{r['ALT']}", axis=1)
konnor_df['ident'] = konnor_df.apply(lambda r: f"{r['POS']}|{r['REF']}|{r['ALT']}", axis=1)


steve_df['QUAL'] = pd.to_numeric(steve_df['QUAL'], 'coerce')
konnor_df['QUAL'] = pd.to_numeric(konnor_df['QUAL'], 'coerce')
steve_df = steve_df[steve_df['QUAL']>30]
konnor_df = konnor_df[konnor_df['QUAL']>30]

setdict = {}
for genome in illumina_samples:
    steve_set = set(steve_df[steve_df['genome_name']==genome]['ident'].to_list())
    konnor_set = set(konnor_df[konnor_df['genome_name']==genome]['ident'].to_list())

    setdict[genome] = steve_set

    genome_intersection = konnor_set.intersection(steve_set)
    konnor_only = konnor_set.difference(steve_set)
    steve_only = steve_set.difference(konnor_set)

    print(f"{genome}")
    print(f"genome_intersection: {len(genome_intersection)}")
    # print(genome_intersection)
    print(f"konnor_only: {len(konnor_only)}")
    # print(konnor_only)
    print(f"steve_only: {len(steve_only)}")
    # print(steve_only)
    # break


control = ["N2_A1A_before_dark_A", "N2_A1A_before_dark_B"]
pheno = ["N2_A1A_afterT1_A",
        "N2_A1A_afterT1_B",
        "N2_A1A_afterT6_A",
        "N2_A1A_afterT6_B",]

control_snps = set()
for g in control:
    control_snps = (control_snps | setdict[g])

pheno_consensus_snps = set.intersection(*[s for g, s in setdict.items() if g in pheno])

just_pheno_snps = pheno_consensus_snps.difference(control_snps)

print(just_pheno_snps)