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
#SBATCH --array=1-23
#
#SBATCH -o logs/metabat2/%j_%a_slurm_output.txt
#SBATCH -e logs/metabat2/%j_%a_slurm_error.txt


wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly"
illumina_adapters="/nfs/chisholmlab001/kve/genomic_resources/references/illumina/all_illumina_adapters.fa"

reads_dir=${wd}/input_data/reads
trimmed_reads_dir=${wd}/intermediate/trimmed_reads

mkdir -p ${trimmed_reads_dir}

read_group=( $(ls ${reads_dir} | grep -o "^[^_]*_[^_]*" | uniq ) )

strain=${read_group[$SLURM_ARRAY_TASK_ID]}
echo $strain


strain_dir=${wd}/"metagenomic_binning"/${strain}
mkdir -p ${strain_dir}

contigs=${wd}/"meta_assemblies"/${strain}/"contigs.fasta"
sam_out=${strain_dir}/"${strain}_mapped.sam"
bam_io=${strain_dir}/"${strain}_mapped.bam"
depth=${strain_dir}/"${strain}_covdepth.txt"
bins_dir=${strain_dir}/"${strain}_bins"
mkdir -p ${bins_dir}

source activate bowtie2

bowtie2-build ${contigs} ${contigs} --threads ${SLURM_NTASKS}
bowtie2 -x ${contigs} -p ${SLURM_NTASKS} -1 ${trimmed_reads_dir}/${strain}_1_sequence.fastq -2 ${trimmed_reads_dir}/${strain}_2_sequence.fastq -S ${sam_out}

conda deactivate

source activate samtools

samtools view -S -b ${sam_out} > ${bam_io}
rm ${sam_out}
samtools sort ${bam_io} -o ${bam_io} 

conda deactivate

source activate metabat2

jgi_summarize_bam_contig_depths --outputDepth ${depth} ${bam_io}

metabat2 -i ${contigs} -a ${depth} -o ${bins_dir}

conda deactivate