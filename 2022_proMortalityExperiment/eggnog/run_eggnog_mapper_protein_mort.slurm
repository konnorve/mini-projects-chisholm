#!/usr/bin/env bash
# 
#SBATCH --job-name=emapPRO
#
# Number of tasks/cores for job
#SBATCH -n 10
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 36:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH --comment="eggNOG mapper"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# Request one node:
#SBATCH --nodes=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate eggnog-mapper

NAME=HTCC7211
IO_DIR=/nfs/chisholmlab001/kve/genomic_resources/strains/het/HTCC7211
IO_GFF=${IO_DIR}/"HTCC7211.gff"
INPUT_FASTA=${IO_DIR}/"HTCC7211_extracted_Proteins.faa"

emapper.py \
--override \
-i ${INPUT_FASTA} \
--itype proteins \
--data_dir /nfs/chisholmlab001/database_repository/eggNOG_Jan_2022 \
-m diamond \
--evalue 0.001 \
--score 60 \
--pident 40 \
--query_cover 20 \
--subject_cover 20 \
--dmnd_ignore_warnings \
--tax_scope auto \
--target_orthologs all \
--report_orthologs \
--go_evidence all \
--pfam_realign none \
--output ${NAME} \
--output_dir ${IO_DIR} \
--temp_dir ${IO_DIR} \
--decorate_gff ${IO_GFF}