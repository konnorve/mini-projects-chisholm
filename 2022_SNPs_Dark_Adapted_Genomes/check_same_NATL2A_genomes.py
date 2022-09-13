from Bio import SeqIO
steve_NATL2A = list(SeqIO.parse("/nfs/chisholmlab001/sbiller/archived_projects/allison_dark_phenotype/NATL2A_genbank.gbk", "genbank"))
konnor_NATL2A = list(SeqIO.parse("/nfs/chisholmlab001/kve/genomic_resources/strains/pro/NATL2A/NATL2A.fna", "fasta"))

print(steve_NATL2A, konnor_NATL2A, sep='\n')
print(steve_NATL2A[0].seq == konnor_NATL2A[0].seq) # evals to true