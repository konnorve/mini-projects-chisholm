#!/usr/bin/env bash
# 
#SBATCH --job-name=bowtie2
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 10-0
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH -n 20
#SBATCH --mem 250G
#
#SBATCH -o logs/%j_%a_slurm_output.txt
#SBATCH -e logs/%j_%a_slurm_error.txt
#
#SBATCH --array=0-65%5

source activate bowtie2

read_dir=/nobackup1/chisholmlab/kve/2022_scope_gradients_HL_adaptation/scratch/1_trimmed_reads

arr=( $(ls ${read_dir}/*_1_trimmed.fastq.gz | sed -E 's/(.*)_1_trimmed.fastq.gz/\1/g' | uniq ))

read_stem=${arr[$SLURM_ARRAY_TASK_ID]}
sample_name=$(basename ${read_stem})

echo $sample_name

#map reads
r1=${read_stem}_1_trimmed.fastq.gz
r2=${read_stem}_2_trimmed.fastq.gz
ru=${read_stem}_unpaired_trimmed.fastq.gz
out_bam=/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/1_mapped_sorted_bams_small_refset/${sample_name}.bam

ls $r1
ls $r2
ls $ru

bowtie2 --no-unal --local -D 15 -R 2 -L 15 -N 1 --gbar 1 --mp 3 -x /nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/small_index/scope_gradients -1 ${r1} -2 ${r2} -U ${ru} -p ${SLURM_NTASKS} | samtools sort -@ ${SLURM_NTASKS} -T ${out_bam}.temp -o ${out_bam}
