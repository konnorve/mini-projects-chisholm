#!/bin/bash
# 
#SBATCH --job-name=busco
#
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=250G
#
# wall clock limit:
#SBATCH --time 24:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm 
#
#SBATCH --comment="busco on rnaSPAdes transcriptome"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH -o logs/busco/%j_slurm_output.txt
#SBATCH -e logs/busco/%j_slurm_error.txt

input="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/transcriptome_qc/cd_hit_est_output/nr_transcriptome.fasta"
output="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/transcriptome_qc/busco_output"

busco_download_path="/nobackup1b/users/chisholmlab/databases_kve/busco"

source activate busco

busco -i ${input} -m "transcriptome" -o ${output} --auto-lineage-euk --download_path ${busco_download_path}
