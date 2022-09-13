import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

gene_dict = pd.read_table("/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_PCC_7942/from_paper/natural_competance_genes.tsv", dtype='str', names=['ID', 'name']).set_index('ID').squeeze("columns").to_dict()

print(gene_dict)

syn_11901_proteins_path = "/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_11901/GCA_005577135.1_ASM557713v1_protein.faa"
syn_8101_proteins_path = "/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_8101/GCF_004209775.1_ASM420977v1_protein.faa"
syn_PCC_7942_proteins_path = "/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_PCC_7942/from_paper/natural_competance_genes.faa"

syn_11901_proteins = {rec.id : rec.seq for rec in SeqIO.parse(syn_11901_proteins_path, "fasta")}
syn_8101_proteins = {rec.id : rec.seq for rec in SeqIO.parse(syn_8101_proteins_path, "fasta")}
syn_PCC_7942_proteins = {rec.id : rec.seq for rec in SeqIO.parse(syn_PCC_7942_proteins_path, "fasta")}

syn_11901_path = "/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_11901/blast_results.tsv"
syn_8101_path = "/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_8101/blast_results.tsv"
syn_PCC_7942_path = "/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_PCC_7942/blast_results.tsv"

def transform_blast_results(df_path):
    blast_table_header = ["query acc.ver",
        "subject acc.ver",
        "% identity",
        "alignment length",
        "mismatches",
        "gap opens",
        "q. start",
        "q. end",
        "s. start",
        "s. end",
        "evalue",
        "bit score",
    ]
    df = pd.read_table(df_path, comment='#', names=blast_table_header)
    df = df[df["% identity"]>0.15]
    df = df.loc[df.groupby("query acc.ver")['bit score'].idxmax()]
    df = df.set_index("query acc.ver")
    return df

syn_11901_blast_df = transform_blast_results(syn_11901_path)
syn_8101_blast_df = transform_blast_results(syn_8101_path)
syn_PCC_7942_blast_df = transform_blast_results(syn_PCC_7942_path)


identities = pd.concat([df['% identity'] for df in [syn_PCC_7942_blast_df, syn_11901_blast_df, syn_8101_blast_df]], axis=1, keys=['syn_PCC_7942', 'syn_11901', 'syn_8101'])
identities = identities / 100
identities = identities.transpose()

fig = plt.figure(figsize=(12, 4), dpi=300)
ax = plt.gca()
plt.imshow(identities, cmap='rainbow')
plt.colorbar()
ax.set_xticks(range(len(identities.columns)), labels=[gene_dict[col.split('_')[1]] for col in identities.columns])
ax.set_yticks(range(len(identities.index)), labels=identities.index)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

ax.set_title("Natural Competance Gene Homology in Syn")
# fig.tight_layout()
plt.savefig("homology.png")