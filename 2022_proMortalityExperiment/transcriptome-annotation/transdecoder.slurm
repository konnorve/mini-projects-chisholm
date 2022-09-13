#!/bin/bash
# 
#SBATCH --job-name=transdecoder
#
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=250G
#
# wall clock limit:
#SBATCH --time 24:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm 
#
#SBATCH --comment="transdecoder on rnaSPAdes transcriptome"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH -o logs/transdecoder/%j_slurm_output.txt
#SBATCH -e logs/transdecoder/%j_slurm_error.txt

input="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/transcriptome_qc/cd_hit_est_output/nr_transcriptome.fasta"
output_dir="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/transcriptome_annotation/transdecoder"

source activate transdecoder

TransDecoder.LongOrfs -t ${input} --output_dir ${output_dir}
TransDecoder.Predict -t ${input} --output_dir ${output_dir}