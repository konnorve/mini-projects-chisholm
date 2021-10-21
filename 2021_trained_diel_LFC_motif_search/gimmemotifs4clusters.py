from gimmemotifs.denovo import gimme_motifs
from pathlib import Path
import sys, time

genome_path = Path('/pool001/kve/2021_trained_diel_LFC_motif_search/NATL2A_genome_references/onlyNATL2A.fna')

proj_dir = Path(sys.argv[1])
num_clusters = int(sys.argv[2])

motif_fastas = proj_dir / 'results' / 'motif_fastas'

gimme_out = proj_dir / 'results' / 'gimme_motifs_results'

n_clust_dir = motif_fastas / '{:02}_clusters'.format(num_clusters)

for motif_fasta in n_clust_dir.iterdir():

    cluster = int(motif_fasta.stem.split('_')[1])

    print(f"cluster {cluster} started from {motif_fasta}")
    start_time = time.time()

    n_cluster_motif_dir = gimme_out / '{:02}_clusters'.format(num_clusters) / 'cluster_{:02}'.format(cluster)
    n_cluster_motif_dir.mkdir(parents=True, exist_ok=True)
    try:
        gimme_motifs(str(motif_fasta), str(n_cluster_motif_dir), params={"genome":genome_path})
    except:
        pass

    td = time.time() - start_time

    print(f"cluster {cluster} finished in {td // 60} minutes")
