#!/usr/bin/env bash
# 
#SBATCH --job-name=trim_reads
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
#SBATCH --comment="trim_reads"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH --mem=250G
#
#SBATCH -o logs/snakemake/%j_slurm_output.txt
#SBATCH -e logs/snakemake/%j_slurm_error.txt

source activate snakemake

snakemake --cores all --resources mem_mb=250000 --use-conda --conda-frontend conda --rerun-incomplete
