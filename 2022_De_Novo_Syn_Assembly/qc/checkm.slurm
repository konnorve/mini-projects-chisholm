#!/usr/bin/env bash
# 
#SBATCH --job-name=syncheckM
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 96:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# Number of tasks/cores for job
#SBATCH -n 1
#
#SBATCH --comment="CheckM on all assemblies"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH --mem 250G
#
#
#SBATCH -o logs/checkm/%j_slurm_output.txt
#SBATCH -e logs/checkm/%j_slurm_error.txt

wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/metagenomic_binning"

source activate checkm

for ASSEMBLY_DIR in ${wd}/*/*_bins/*
do
    echo $ASSEMBLY_DIR
    CHECKM_DIR=${ASSEMBLY_DIR}/checkM
    checkm lineage_wf -t ${SLURM_NTASKS} -x fa ${ASSEMBLY_DIR} ${CHECKM_DIR}
done