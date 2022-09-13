import pandas as pd
import os
from pathlib import Path
import shutil

wd = Path('/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database')
df = pd.read_table(wd / "lists" / "old" / "MARMICRODB_catalog.tsv")
genomes_dir = Path(wd / 'marmicro_hets_cyanos_no_prosyn' / 'rest')

df = df.join(df['lineage_assignment'].str.split(';', expand=True))

df = df[df[1]=='Proteobacteria']

for species in df[2].unique():
    species = species.replace("/", "_").replace(" ", "_")
    print(species)
    subset = df[df[2]==species]
    species_dir = wd / "marmicro_hets_cyanos_no_prosyn" / species
    species_genome_dir = species_dir / "refs"
    species_fastani_dir = species_dir / "fastani_outputs"
    species_dir.mkdir(exist_ok=True)
    species_genome_dir.mkdir(exist_ok=True)
    species_fastani_dir.mkdir(exist_ok=True)
    genomes = [f"{Path(f).name}_genomic.fna.gz" for f in subset.assembly_ftp.unique() if not isinstance(f, float)]
    for genome in genomes:
        try:
            shutil.move(genomes_dir / genome, species_genome_dir / genome)
        except FileNotFoundError as e:
            pass
    full_paths = [str(f) for f in species_genome_dir.iterdir()]
    with open(species_dir / f"{species}.txt", "w") as f:
        f.write("\n".join(map(str, full_paths)))
    subset.to_csv(species_dir / f"{species}.tsv", sep='\t', index=False)
