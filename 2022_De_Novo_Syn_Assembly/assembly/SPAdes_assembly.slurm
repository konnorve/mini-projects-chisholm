#!/usr/bin/env bash
# 
#SBATCH --job-name=SPAdes
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
#SBATCH --array=0-23
#
#SBATCH -o logs/SPAdes/%j_%a_slurm_output.txt
#SBATCH -e logs/SPAdes/%j_%a_slurm_error.txt

wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly"
illumina_adapters="/nfs/chisholmlab001/kve/genomic_resources/references/illumina/all_illumina_adapters.fa"

reads_dir=${wd}/input_data/reads
trimmed_reads_dir=${wd}/intermediate/trimmed_reads

mkdir -p ${trimmed_reads_dir}

read_group=( $(ls ${reads_dir} | grep -o "^[^_]*_[^_]*" | uniq ) )

strain=${read_group[$SLURM_ARRAY_TASK_ID]}
echo $strain

source activate bbmap

bbduk.sh \
    in1=${reads_dir}/${strain}_1_sequence.fastq in2=${reads_dir}/${strain}_2_sequence.fastq \
    out1=${trimmed_reads_dir}/${strain}_1_sequence.fastq out2=${trimmed_reads_dir}/${strain}_2_sequence.fastq \
    minlen=25 qtrim=rl trimq=10 \
    ref=${illumina_adapters} ktrim=r k=23 mink=11 hdist=1

conda deactivate

uniq_assem_dir=${wd}/"assemblies"/${strain}
mkdir -p ${uniq_assem_dir}
/home/kve/bin/SPAdes-3.15.4-Linux/bin/spades.py --isolate -t ${SLURM_NTASKS} -1 ${trimmed_reads_dir}/${strain}_1_sequence.fastq -2 ${trimmed_reads_dir}/${strain}_2_sequence.fastq -o ${uniq_assem_dir}

uniq_assem_dir=${wd}/"meta_assemblies"/${strain}
mkdir -p ${uniq_assem_dir}
/home/kve/bin/SPAdes-3.15.4-Linux/bin/spades.py --meta -t ${SLURM_NTASKS} -1 ${trimmed_reads_dir}/${strain}_1_sequence.fastq -2 ${trimmed_reads_dir}/${strain}_2_sequence.fastq -o ${uniq_assem_dir}




