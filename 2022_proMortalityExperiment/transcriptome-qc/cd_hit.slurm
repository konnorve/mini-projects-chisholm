#!/bin/bash
# 
#SBATCH --job-name=cdhit
#
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=250G
#
# wall clock limit:
#SBATCH --time 24:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm 
#
#SBATCH --comment="cdhit on rnaSPAdes transcriptome"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH -o logs/transrate/%j_slurm_output.txt
#SBATCH -e logs/transrate/%j_slurm_error.txt

transcriptome_dir="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/rnaSPAdes"

transcriptome=${transcriptome_dir}/"transcripts.fasta"

read_dir="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/rnaSPAdes/concatenated_reads/"

left_reads=${read_dir}/"all_1_trimmed.fastq.gz"
right_reads=${read_dir}/"all_2_trimmed.fastq.gz"

output_file="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/transcriptome_qc/cd_hit_est_output/nr_transcriptome.fasta"

source activate cd-hit

cd-hit-est -i ${transcriptome} -M 200000 -n 8 -c 0.9 -o ${output_file} -T ${SLURM_NTASKS}

