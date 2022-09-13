#!/usr/bin/env bash
# 
#SBATCH --job-name=mergealn
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

source activate phykit

wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/extract_align_gtdbtk_bac120/msa_bac120"

cd ${wd}

for f in gene_msas/*.aln
do
    realpath $f >> alignment_files.txt
done

phykit create_concat -a alignment_files.txt -p all_syn_bac120_msa

tr ? -  < all_syn_bac120_msa.fa > all_syn_bac120_msa_for_tree.fa

awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' all_syn_bac120_msa_for_tree.fa