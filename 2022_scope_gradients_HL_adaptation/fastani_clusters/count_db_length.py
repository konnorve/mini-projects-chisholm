from Bio.SeqIO.FastaIO import SimpleFastaParser
from pathlib import Path

total_len = 0
for i, genome in enumerate(Path('/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/chosen_references').iterdir()):
    print(i, genome, sep='\t')
    if genome.suffix == '.fna':
        with open(genome) as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                total_len += len(seq)
print(total_len)