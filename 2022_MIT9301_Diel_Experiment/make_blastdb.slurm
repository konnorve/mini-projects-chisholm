#!/usr/bin/env bash
# 
#SBATCH --job-name=blastdb
#
# Number of tasks/cores for job
#SBATCH -n 20
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 24:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH --comment="make blastdb"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH --mem=250G
#
#SBATCH -o logs/blastdb/%j_slurm_output.txt
#SBATCH -e logs/blastdb/%j_slurm_error.txt

source activate blast

for file in /nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/comparison_seqs/*
do
    makeblastdb -in ${file}
done

