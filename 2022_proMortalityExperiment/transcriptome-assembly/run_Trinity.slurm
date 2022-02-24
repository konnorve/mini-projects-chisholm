#!/usr/bin/env bash
# 
#SBATCH --job-name=trinity
#
# Number of tasks/cores for job
#SBATCH -n 10
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
#SBATCH --comment="eggNOG mapper"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# Request one node:
#SBATCH --nodes=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
#SBATCH --mem=250G
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate trinity_helpers

SAMPLE_FILE=/home/kve/scripts/mini_projects/2022_proMortalityExperiment/transcriptome-assembly/polyAAA_sample_table.tsv

OUT_DIR=/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/trinity_output

/home/kve/bin/trinityrnaseq-v2.13.2/Trinity \
    --trimmomatic \
    --seqType fq --max_memory 250G --CPU 10 \
    --samples_file ${SAMPLE_FILE} \
    --output ${OUT_DIR}
    