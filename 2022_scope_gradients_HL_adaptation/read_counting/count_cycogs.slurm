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
#SBATCH --array=0

readarray -t arr < 'bams.txt'
htseq_out_dir=${arr[$SLURM_ARRAY_TASK_ID]}
sample_name=$(basename ${htseq_out_dir})
bam_stem="/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/1_mapped_sorted_bams/${sample_name}"
bam="${bam_stem}.bam"
echo $sample_name

mkdir -p $htseq_out_dir

#map reads
bam_se=${bam_stem}_unpaired.bam
bam_pe=${bam_stem}_paired.bam

source activate samtools
if [[ ! -f ${bam_se} ]]; then samtools view -F1 -b -o ${bam_se} ${bam}; fi
if [[ ! -f ${bam_pe} ]]; then samtools view -f1 -b -o ${bam_pe} ${bam}; fi

conda deactivate
source activate htseq

for gff in /nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/chosen_reference_gffs/*.gff
do
    genome=$(basename ${gff})
    genome=${genome%.*}
    output=${htseq_out_dir}/${genome}.tsv
    if [[ ! -f ${output} ]]
    then
        echo $genome
        htseq-count --idattr=gene -a 30 -s no -t CDS -r pos --nonunique all -c ${output} ${bam_pe} ${bam_se} ${gff}
    fi
done

rm ${bam_pe} ${bam_se}