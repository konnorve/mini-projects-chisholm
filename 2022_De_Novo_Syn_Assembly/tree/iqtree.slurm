#!/usr/bin/env bash
# 
#SBATCH --job-name=iqtree
#
# wall clock limit:
#SBATCH --time 10-0
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH --nodes=1
#SBATCH -n 20
#SBATCH --mem 250G
#
#SBATCH -o ../logs/iqtree/%j_%a_slurm_output.txt
#SBATCH -e ../logs/iqtree/%j_%a_slurm_error.txt


MSA="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/extract_align_gtdbtk_bac120/msa_bac120/all_syn_bac120_msa_for_tree.fa"
OUTDIR="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/iqtree"

/home/kve/bin/iqtree-1.6.12-Linux/bin/iqtree -s ${MSA} -m "LG+F+R10" -nt ${SLURM_NTASKS} -pre ${OUTDIR}/all_syn_iqtree

