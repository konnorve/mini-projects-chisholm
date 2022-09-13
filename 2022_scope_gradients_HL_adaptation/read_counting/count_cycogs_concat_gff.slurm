#!/usr/bin/env bash
# 
#SBATCH --job-name=htseq
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
#SBATCH -n 1
#SBATCH --mem 10G
#
#SBATCH -o logs/%j_%a_slurm_output.txt
#SBATCH -e logs/%j_%a_slurm_error.txt
#
#SBATCH --array=0-65

gff="/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/concat_reference_annotations.gff"

readarray -t arr < 'bams.txt'
bam=${arr[$SLURM_ARRAY_TASK_ID]}
bam_stem=${bam%.*}
sample_name=$(basename ${bam_stem})
htseq_out="/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/2_htseq_concat_output_0_mapq/${sample_name}.tsv"
echo $sample_name

#map reads
bam_se=${bam_stem}_unpaired.bam
bam_pe=${bam_stem}_paired.bam

source activate samtools
if [[ ! -f ${bam_se} ]]; then samtools view -F1 -b -o ${bam_se} ${bam}; fi
if [[ ! -f ${bam_pe} ]]; then samtools view -f1 -b -o ${bam_pe} ${bam}; fi

conda deactivate
source activate htseq

htseq-count --idattr=gene -a 0 -s no -t CDS -r pos --nonunique all -c ${htseq_out} ${bam_pe} ${bam_se} ${gff}

# rm ${bam_pe} ${bam_se}