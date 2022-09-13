import pandas as pd
from igraph import Graph
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pathlib import Path
import gzip
import yaml
import shutil

wd = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/marmicro_hets_cyanos_no_prosyn/Alphaproteobacteria/")

chosen_references_dir = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/chosen_references")

df = pd.read_table(wd / "all_fastani_results.txt", names=['q','r','ani','fragment_mappings','query_fragments'])

print(df)

g = Graph()

unique_genomes = list(df.q.unique())

g.add_vertices(len(unique_genomes))
g.vs['name'] = unique_genomes

for i, q in enumerate(unique_genomes):
    for r in df[(df.q==q) & (df.ani>95)].r.values:
        g.add_edges([(i, unique_genomes.index(r))])

vc = g.community_leiden()

clusters = {}
for i, cluster in enumerate(iter(vc)):
    print(i, cluster)
    seq_lengths = {}
    for j in cluster:
        genome = unique_genomes[j]
        print(genome)
        with gzip.open(genome, 'rt') as in_handle:
            total_len = 0
            for title, seq in SimpleFastaParser(in_handle):
                total_len += len(seq)
            seq_lengths[genome] = total_len
    representative_genome = Path(max(seq_lengths, key=seq_lengths.get))
    shutil.copy(representative_genome, chosen_references_dir)
    clusters[i] = seq_lengths

with open(wd / 'clusters.txt', 'w') as clusters_handle:
    yaml.dump(clusters, clusters_handle)