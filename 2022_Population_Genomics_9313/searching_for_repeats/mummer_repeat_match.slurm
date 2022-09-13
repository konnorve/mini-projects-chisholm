#!/usr/bin/env bash
# 
#SBATCH --job-name=repeat
#
# wall clock limit:
#SBATCH --time 1-0
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem 250000
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate mummer

potential_duplication="/nfs/chisholmlab001/kve/2022_Population_Genomics_9313/genome_duplication_analysis/potential_duplication.fasta"
repeat-match -E ${potential_duplication} > /nfs/chisholmlab001/kve/2022_Population_Genomics_9313/genome_duplication_analysis/potential_duplication.repeats

original_genome="/nfs/chisholmlab001/kve/genomic_resources/strains/pro/MIT9313/MIT9313.fna"
repeat-match -E ${original_genome} > /nfs/chisholmlab001/kve/2022_Population_Genomics_9313/genome_duplication_analysis/original_genome.repeats
