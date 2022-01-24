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

# NAME=HTCC7211_eggNOG_Prodigal
# OUT_DIR=/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations/HTCC7211_Genomic_Resources/${NAME}
# INPUT_FASTA=/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/input_data/culture_genome_refs/HTCC7211.fna

NAME=MIT9301_eggNOG_Prodigal
OUT_DIR=/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations/MIT9301_Genomic_Resources/${NAME}
INPUT_FASTA=/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/input_data/culture_genome_refs/MIT9301.fna

mkdir -p ${OUT_DIR}

emapper.py \
--data_dir /nfs/chisholmlab001/database_repository \
--output ${NAME} \
--output_dir ${OUT_DIR} \
--temp_dir ${OUT_DIR} \
--override \
-m diamond \
--dmnd_ignore_warnings \
-i ${INPUT_FASTA} \
--evalue 0.001 \
--score 60 \
--pident 40 \
--query_cover 20 \
--subject_cover 20 \
--itype genome \
--genepred prodigal \
--tax_scope auto \
--target_orthologs all \
--go_evidence non-electronic \
--pfam_realign none \
--report_orthologs \
--decorate_gff yes # > ${OUT_DIR}/emapper.out 2> ${OUT_DIR}/emapper.err