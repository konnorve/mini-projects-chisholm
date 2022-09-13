#!/usr/bin/env bash
# 
#SBATCH --job-name=bt2idx
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 10-0
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH -n 20
#SBATCH --mem 250G
#
#SBATCH -o logs/%j_%a_slurm_output.txt
#SBATCH -e logs/%j_%a_slurm_error.txt

source activate bowtie2

cd /nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database

file_list=$(ls chosen_references_small/*.fna | tr '\n' ',')

echo $file_list
echo

bowtie2-build --threads ${SLURM_NTASKS} $file_list small_index/scope_gradients