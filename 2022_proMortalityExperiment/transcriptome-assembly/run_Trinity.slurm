#!/bin/bash
# 
#SBATCH --job-name=trinity
#
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=500G
#
# wall clock limit:
#SBATCH --time 12:00:00
#
# Partition to use:
#SBATCH --partition newnodes 
#
#SBATCH --comment="trinity transcriptome assembly of grazer"
#
#SBATCH --constraint=centos7
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate trinity-dependencies

module load trinityrnaseq/2.14.0

SAMPLE_FILE=/home/kve/scripts/mini_projects/2022_proMortalityExperiment/transcriptome-assembly/polyAAA_sample_table.tsv

OUT_DIR=/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/trinity_output

Trinity \
    --trimmomatic \
    --seqType fq \
    --bflyHeapSpaceInit 10G \
    --max_memory 240G \
    --CPU ${SLURM_NTASKS} \
    --samples_file ${SAMPLE_FILE} \
    --output ${OUT_DIR}
    