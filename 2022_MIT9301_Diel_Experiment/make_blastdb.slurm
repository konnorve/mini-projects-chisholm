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

cd /nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/MIT9301_seq_db/

makeblastdb -in "MIT9301_IMG_2623620961.fna" -dbtype nucl -max_file_sz 4GB -out "MIT9301_IMG_2623620961_seq"
