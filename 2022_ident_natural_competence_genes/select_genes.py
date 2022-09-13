import pandas as pd
from Bio import SeqIO

gene_subset = pd.read_table("/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_PCC_7942/from_paper/natural_competance_genes.tsv", names=['ID', 'name'])

all_proteins = list(SeqIO.parse("/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_PCC_7942/original_genome/genes.faa", "fasta"))

natural_competance_ids = [f"Synpcc7942_{i:04d}" for i in gene_subset['ID']]

natural_competance_gene_list = [p for p in all_proteins if p.id in natural_competance_ids]

print(len(natural_competance_gene_list))
print(natural_competance_gene_list)

SeqIO.write(natural_competance_gene_list, "/nfs/chisholmlab001/kve/2022_ident_natural_competence_genes/genomic_resources/syn_PCC_7942/from_paper/natural_competance_genes.faa", "fasta")