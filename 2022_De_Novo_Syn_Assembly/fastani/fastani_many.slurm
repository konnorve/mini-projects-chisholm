#!/usr/bin/env bash
# 
#SBATCH --job-name=fastani
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 96:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# Number of tasks/cores for job
#SBATCH -n 20
#
#SBATCH --comment="fastani on de novo syn isolates"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# Request nodes:
#SBATCH --nodes=1
#
#SBATCH --mem 250G
#
#SBATCH -o ../logs/fastani/%j_slurm_output.txt
#SBATCH -e ../logs/fastani/%j_slurm_error.txt

source activate fastani

fastANI --ql syn_genomes.txt --rl syn_genomes.txt -o fastANI.txt -t ${SLURM_NTASKS} --matrix