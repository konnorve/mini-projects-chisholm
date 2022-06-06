#!/usr/bin/env bash
# 
#SBATCH --job-name=phylosift
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
#SBATCH --comment="phylosift on de novo syn isolates"
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
#SBATCH -o ../logs/phylosift/%j_slurm_output.txt
#SBATCH -e ../logs/phylosift/%j_slurm_error.txt

wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly"

arr=(${wd}/all_syn_fastas/*)

fasta_input=${arr[0]}
output_dir=$(basename $fasta_input)
output_dir=${wd}/phylosift_results/${output_dir%.fa}

mkdir -p $output_dir

/home/kve/bin/phylosift_v1.0.1/bin/phylosift all --threads ${SLURM_NTASKS} --force --output ${output_dir} ${fasta_input}