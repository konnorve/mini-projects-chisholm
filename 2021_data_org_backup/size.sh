#!/usr/bin/env bash
# 
#SBATCH --job-name=size_check
#
# Number of tasks/cores for job
#SBATCH -n 1
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 36:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH --comment="Check size of key directories"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
#SBATCH -o %j_slurm_output.txt
#SBATCH -e %j_slurm_error.txt

for DATA_DIR in /nobackup1/kve # /nobackup1/thomase2 /nobackup1/thackl /nobackup1/chisholmlab
do
DIR_BASENAME=$(basename $DATA_DIR)
du -h -d 2 $DATA_DIR >> ~/data_backup_logistics/data_amount_$DIR_BASENAME.txt
done