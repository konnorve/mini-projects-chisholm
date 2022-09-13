#!/usr/bin/env bash
# 
#SBATCH --job-name=fastani
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 60
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm,newnodes,sched_mit_hill,sched_any
#
#SBATCH -n 10
#SBATCH --mem 120G
#
#SBATCH -o logs/%j_%a_slurm_output.txt
#SBATCH -e logs/%j_%a_slurm_error.txt
#
#SBATCH --array=1500-1845 # 1845

source activate fastani

wd="/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/marmicro_hets_cyanos_no_prosyn"
wd=${wd}/"Alphaproteobacteria"
genome_list=${wd}/"Alphaproteobacteria.txt"

mapfile -t arr < ${genome_list}

q_genome=${arr[$SLURM_ARRAY_TASK_ID]}

genome_name=$(basename $q_genome)
genome_name=${genome_name%.fna.gz}
output_file=${wd}/"fastani_outputs"/${genome_name}.fastani.txt

echo $q_genome
echo $genome_name
echo $output_file

fastANI -q ${q_genome} --rl ${genome_list} --output ${output_file} -t ${SLURM_NTASKS}