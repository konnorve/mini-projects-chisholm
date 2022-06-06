#!/bin/bash
# 
#SBATCH --job-name=transRate
#
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=250G
#
# wall clock limit:
#SBATCH --time 12:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm 
#
#SBATCH --comment="transRate on rnaSPAdes transcriptome"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH -o logs/transrate/%j_slurm_output.txt
#SBATCH -e logs/transrate/%j_slurm_error.txt

transcriptome_dir="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/rnaSPAdes"

hard_transcriptome=${transcriptome_dir}/"hard_filtered_transcripts.fasta"
firm_transcriptome=${transcriptome_dir}/"transcripts.fasta"
soft_transcriptome=${transcriptome_dir}/"soft_filtered_transcripts.fasta"

read_dir="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/rnaSPAdes/concatenated_reads/"

left_reads=${read_dir}/"all_1_trimmed.fastq.gz"
right_reads=${read_dir}/"all_2_trimmed.fastq.gz"

output_dir="/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/transcriptome_qc/transrate_og_rnaSPAdes_assembly"

source activate transrate

transrate --assembly ${hard_transcriptome},${firm_transcriptome},${soft_transcriptome} --left ${left_reads} --right ${right_reads} --threads ${SLURM_NTASKS} --output ${output_dir}

