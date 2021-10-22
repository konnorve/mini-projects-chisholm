#!/usr/bin/env bash
# 
#SBATCH --job-name=LFC_clust
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 48:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# Number of tasks/cores for job
#SBATCH -n 10
#
#SBATCH --comment="LFC Clustering analysis"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# Request nodes:
#SBATCH --nodes=1
#
#SBATCH --mem 250000
#
#SBATCH -o ../logs/%j_%a_slurm_output.txt
#SBATCH -e ../logs/%j_%a_slurm_error.txt

source activate 2021_trained_diel_LFC_motif_search

python3 determine_dark_tolerant_motifs.py "/pool001/kve/2021_trained_diel_LFC_motif_search"