#!/usr/bin/env bash
# 
#SBATCH --job-name=gtdb
#
# wall clock limit:
#SBATCH --time 96:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH -o logs/gtdbtk/%j_slurm_output.txt
#SBATCH -e logs/gtdbtk/%j_slurm_error.txt

cd /nobackup1/chisholmlab/databases_kve

wget -nv https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz

tar -xvzf gtdbtk_v2_data.tar.gz

