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
#SBATCH --array=0-23
#
#SBATCH -o ../logs/gtdbtk/%j_%a_slurm_output.txt
#SBATCH -e ../logs/gtdbtk/%j_%a_slurm_error.txt

wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/metagenomic_binning"
read_group=( $(ls ${wd} ) )
strain=${read_group[$SLURM_ARRAY_TASK_ID]}
echo $strain

genome_dir=${wd}/${strain}/"bin_fastas"
out_dir=${wd}/${strain}/"gtdbtk_identify"

source activate gtdbtk

gtdbtk identify --genome_dir ${genome_dir} --extension "fa" --out_dir ${out_dir} --prefix ${strain} --cpus ${SLURM_NTASKS} --write_single_copy_genes