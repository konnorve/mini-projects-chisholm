#!/usr/bin/env bash
# 
#SBATCH --job-name=bwa
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 4-0
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --mem 250G
#
#
#SBATCH -o logs/bwa/%j_slurm_output.txt
#SBATCH -e logs/bwa/%j_slurm_error.txt

source activate bwa

wd="/nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/artifical_read_mapping"

ref_genome=${wd}/joint_MIT9301_internal_standards.fna
read_1=${wd}/"hiseq2000_art_reads"/"art_MIT9301_1.fq"
read_2=${wd}/"hiseq2000_art_reads"/"art_MIT9301_2.fq"
sam=${wd}/"hiseq2000_mapping"/"mapped_MIT9301.sam"
bam_unsort=${wd}/"hiseq2000_mapping"/"mapped_MIT9301_unsorted.bam"
bam_sort=${wd}/"hiseq2000_mapping"/"mapped_MIT9301_sorted.bam"
contig_info=${wd}/"hiseq2000_mapping"/"MIT9301_mapping_info.txt"

bwa index ${ref_genome}
bwa mem -t ${SLURM_NTASKS} ${ref_genome} ${read_1} ${read_2} > ${sam}

conda deactivate
source activate samtools

samtools view -@ ${SLURM_NTASKS} -S -b ${sam} > ${bam_unsort}
samtools sort -@ ${SLURM_NTASKS} ${bam_unsort} -o ${bam_sort}
samtools index -b ${bam_sort}
samtools idxstats ${bam_sort} > ${contig_info}

rm ${sam} ${bam_unsort}