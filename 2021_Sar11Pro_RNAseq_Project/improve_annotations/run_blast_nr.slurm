#!/usr/bin/env bash
# 
#SBATCH --job-name=blastp
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 48:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# Number of tasks/cores for job
#SBATCH -n 10
#
#SBATCH --comment="Blast for improving MIT9301 annotations"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# Request nodes:
#SBATCH --nodes=1
#
#SBATCH --mem 250000
#
#SBATCH -o logs/%j_%a_slurm_output.txt
#SBATCH -e logs/%j_%a_slurm_error.txt

source activate blast

MIT9301_CDS="/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations/MIT9301_extracted_cds.faa"
BLAST_OUTPUT="/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations/MIT9301_cds_blastp_nr.xml"

blastp -db nr -remote -query $MIT9301_CDS -outfmt 5 -out $BLAST_OUTPUT

