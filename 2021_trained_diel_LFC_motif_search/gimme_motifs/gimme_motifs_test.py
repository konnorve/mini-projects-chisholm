from gimmemotifs.denovo import gimme_motifs, create_background
from gimmemotifs.scanner import Scanner as s


from pathlib import Path
input_fasta = Path('/pool001/kve/2021_trained_diel_LFC_motif_search/results/motif_fastas/09_clusters/cluster_08.fasta')
output_dir = Path('/pool001/kve/2021_trained_diel_LFC_motif_search/results/gimme_motifs_results/09_clusters/cluster_08')
genome_path = Path('/pool001/kve/2021_trained_diel_LFC_motif_search/NATL2A_genome_references/onlyNATL2A.fna')


params={
    "genome":genome_path
}

motifs = gimme_motifs(str(input_fasta), str(output_dir), params=params)
