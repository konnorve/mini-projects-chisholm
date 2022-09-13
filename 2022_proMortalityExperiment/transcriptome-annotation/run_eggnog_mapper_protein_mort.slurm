#!/usr/bin/env bash
# 
#SBATCH --job-name=eggNOGm
#
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
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=250G
#
#SBATCH -o logs/eggnog/%j_slurm_output.txt
#SBATCH -e logs/eggnog/%j_slurm_error.txt

source activate eggnog-mapper

NAME=nr_transcriptome
IO_DIR=/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/transcriptome_annotation/eggNOGmapper
IO_GFF=${IO_DIR}/"nr_transcriptome.fasta.transdecoder.CDS.gff3"
INPUT_FASTA=${IO_DIR}/"nr_transcriptome_cds.faa"

emapper.py \
--override \
-i ${INPUT_FASTA} \
--itype proteins \
--data_dir /nobackup1b/users/chisholmlab/databases_kve/eggNOG_Jan_2022 \
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
