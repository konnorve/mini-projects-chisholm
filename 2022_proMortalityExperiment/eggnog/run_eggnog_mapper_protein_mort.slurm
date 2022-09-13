#!/usr/bin/env bash
# 
#SBATCH --job-name=emapPRO
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 12:00:00
#
# Partition to use:
#SBATCH --partition newnodes,sched_mit_chisholm
#
#SBATCH --mem 250G
#SBATCH -n 20
#
#SBATCH --comment="eggNOG mapper"
#
#SBATCH --array=0-23
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt


wd="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/input_data_ribozero/het_bins"

bins=( $(ls ${wd}) )

NAME=${bins[$SLURM_ARRAY_TASK_ID]}
echo $NAME

IO_DIR=/nfs/chisholmlab001/kve/2022_proMortalityExperiment/input_data_ribozero/het_bins/$NAME
# IO_GFF=${IO_DIR}/"$NAME.gff"
INPUT_FASTA=${IO_DIR}/"$NAME.fna"


source activate eggnog-mapper

emapper.py \
--override \
-i ${INPUT_FASTA} \
--itype genome \
--genepred prodigal \
--data_dir /nobackup1/chisholmlab/databases_kve/eggNOG_Jan_2022 \
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
--decorate_gff yes \
--cpu ${SLURM_NTASKS}