#!/usr/bin/env bash
# 
#SBATCH --job-name=vibrant
#
# wall clock limit:
#SBATCH --time 5-0
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 50000
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate vibrant

input="/nfs/chisholmlab001/kve/2021_HIFI_Genome_Closing_Project/zzz_pulling_contigs_for_sean_july_11_2022/pacbio_WH7803x.fasta"
outdir="/nfs/chisholmlab001/kve/2022_WH7803_variant_calling/vibrant_WH7803x"
databases="/nobackup1/chisholmlab/databases_kve/vibrantv1.2.1"

VIBRANT_run.py -i ${input} -folder ${outdir} -d ${databases}

