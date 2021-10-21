#!/usr/bin/env bash
# 
#SBATCH --job-name=GimmeMotifs
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
#SBATCH --comment="GimmeMotifs"
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
# 9, 18, 39, 59
#SBATCH --array=9,18,39,59
#
#SBATCH -o logs/%j_%a_slurm_output.txt
#SBATCH -e logs/%j_%a_slurm_error.txt

source activate gimme

python3 gimmemotifs4clusters.py "/pool001/kve/2021_trained_diel_LFC_motif_search" ${SLURM_ARRAY_TASK_ID}