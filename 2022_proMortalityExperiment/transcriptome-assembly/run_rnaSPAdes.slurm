#!/bin/bash
# 
#SBATCH --job-name=rnaSPAdes
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
#SBATCH --comment="rnaSPAdes transcriptome assembly of grazer"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

# SAMPLE_FILE=/home/kve/scripts/mini_projects/2022_proMortalityExperiment/transcriptome-assembly/polyAAA_samples.yaml

SAMPLE_DIR="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/rnaSPAdes/concatenated_reads"

OUT_DIR=/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/rnaSPAdes

# /home/kve/bin/SPAdes-3.12.0-Linux/bin/spades.py --rna -1 ${SAMPLE_DIR}/"all_1_trimmed.fastq.gz" -2 ${SAMPLE_DIR}/"all_2_trimmed.fastq.gz" -t ${SLURM_NTASKS} --cov-cutoff 300 --memory 500 -o ${OUT_DIR}

/home/kve/bin/SPAdes-3.12.0-Linux/bin/spades.py -o ${OUT_DIR} --continue


