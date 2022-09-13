#!/usr/bin/env bash
# 
#SBATCH --job-name=gtdbtk
#
# wall clock limit:
#SBATCH --time 96:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# Request nodes:
#SBATCH --nodes=1
#
# Number of tasks/cores for job
#SBATCH -n 20
#
#SBATCH --mem 250G
#
#SBATCH -o ../logs/gtdbtk/%j_%a_slurm_output.txt
#SBATCH -e ../logs/gtdbtk/%j_%a_slurm_error.txt

wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/"

genome_dir=${wd}/"all_syn_fastas"
out_dir=${wd}/"extract_align_gtdbtk_bac120"/"gtdbtk_identify_bac120"

source activate gtdbtk

gtdbtk identify --genome_dir ${genome_dir} --extension "fna" --out_dir ${out_dir} --cpus ${SLURM_NTASKS} --write_single_copy_genes