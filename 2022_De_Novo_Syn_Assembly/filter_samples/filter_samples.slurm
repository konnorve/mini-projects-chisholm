#!/usr/bin/env bash
# 
#SBATCH --job-name=filter
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
#SBATCH --comment="read filtering using de novo syn isolate bins"
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
#SBATCH --array=0-25
#
#SBATCH -o ../logs/filter_samples/%j_%a_slurm_output.txt
#SBATCH -e ../logs/filter_samples/%j_%a_slurm_error.txt


wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly"

bin_dir=${wd}/all_de_novo_bins
outdir=${wd}/filtered_readsets
sample_dir=${wd}/input_data/reads

bins=( $(ls ${bin_dir}/*.fa) )

bin=${bins[$SLURM_ARRAY_TASK_ID]}
bin=$(basename $bin)
echo $bin

sample=$(echo $bin | sed -r "s/(\S+?)_bins.*/\1/")

bin_contigs=${bin_dir}/${bin}
fwd_reads=${sample_dir}/"${sample}_1_sequence.fastq"
rev_reads=${sample_dir}/"${sample}_2_sequence.fastq"
filtered_fwd_reads=${outdir}/"${bin%.*}_1_sequence.fastq"
filtered_rev_reads=${outdir}/"${bin%.*}_2_sequence.fastq"
sam_out=${outdir}/"${bin}_subreadset_mapped.sam"

source activate bowtie2

bowtie2-build ${bin_contigs} ${bin_contigs} --threads ${SLURM_NTASKS}
bowtie2 --local -x ${bin_contigs} -p ${SLURM_NTASKS} -1 ${fwd_reads} -2 ${rev_reads} -S ${sam_out}

conda deactivate
source activate samtools

samtools fastq -F 4 -@ ${SLURM_NTASKS} ${sam_out} -1 ${filtered_fwd_reads} -2 ${filtered_rev_reads} -0 /dev/null -s /dev/null

rm ${sam_out}