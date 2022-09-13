#!/usr/bin/env bash
# 
#SBATCH --job-name=mafft
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
#SBATCH -o ../logs/mafft/%j_%a_slurm_output.txt
#SBATCH -e ../logs/mafft/%j_%a_slurm_error.txt

bac120_dir="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/extract_align_gtdbtk_bac120/gtdbtk_identify_bac120/identify/intermediate_results/single_copy_fasta/bac120"

gene_msas_dir="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/extract_align_gtdbtk_bac120/msa_bac120/gene_msas"

for fasta in ${bac120_dir}/*
do
    /home/kve/bin/mafft-7.490-with-extensions/bin/mafft --auto ${fasta} > ${gene_msas_dir}/$(basename ${fasta%.*}).aln
done

