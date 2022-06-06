#!/usr/bin/env bash
# 
#SBATCH --job-name=kaiju
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
#SBATCH -n 20
#
#SBATCH --comment="SPAdes assembly de novo syn isolates"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# Request nodes:
#SBATCH --nodes=1
#
#SBATCH --mem 250G
#
#SBATCH -o logs/kaiju/%j_%a_slurm_output.txt
#SBATCH -e logs/kaiju/%j_%a_slurm_error.txt

wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/metagenomic_binning"

kaiju_index_file="/nobackup1b/users/chisholmlab/databases_kve/kaiju-20220425/refseq/kaiju_db_refseq.fmi"
nodes_dmp="/nobackup1b/users/chisholmlab/databases_kve/kaiju-20220425/nodes.dmp"
names_dmp="/nobackup1b/users/chisholmlab/databases_kve/kaiju-20220425/names.dmp"

for assembly_path in $wd/*/*_bins/*/*.fa
do
    assembly_file_name=$(basename "$assembly_path")
    assembly_name="${assembly_file_name%.*}"
    ASSEMBLY_DIR=$(dirname "$assembly_path")

    kaiju_bc_result_txt=${ASSEMBLY_DIR}/"kaiju_classification_output.${assembly_name}.txt"
    kaiju_bc_result_table=${ASSEMBLY_DIR}/"kaiju_classification_output.${assembly_name}.tsv"
    kaiju_bc_result_names=${ASSEMBLY_DIR}/"kaiju_classification_output_named.${assembly_name}.tsv"
    kaiju_bc_result_name_class=${ASSEMBLY_DIR}/"kaiju_classification_output_named_classified.${assembly_name}.tsv"

    echo ${assembly_path}

    /home/kve/bin/kaiju/bin/kaiju -v -t ${nodes_dmp} -f ${kaiju_index_file} -i ${assembly_path} -o ${kaiju_bc_result_txt} -z ${SLURM_NTASKS}

    /home/kve/bin/kaiju/bin/kaiju2table -t ${nodes_dmp} -n ${names_dmp} -r genus -o ${kaiju_bc_result_table} ${kaiju_bc_result_txt}

    /home/kve/bin/kaiju/bin/kaiju-addTaxonNames -t ${nodes_dmp} -n ${names_dmp} -i ${kaiju_bc_result_txt} -o ${kaiju_bc_result_names}

    # removed unclassified contigs
    sed '/^U/d' ${kaiju_bc_result_names} > ${kaiju_bc_result_name_class}

    rm ${kaiju_bc_result_table} ${kaiju_bc_result_txt} ${kaiju_bc_result_names}
done