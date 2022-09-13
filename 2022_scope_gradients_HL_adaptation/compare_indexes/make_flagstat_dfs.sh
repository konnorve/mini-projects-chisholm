source activate samtools

for file in /nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/1_mapped_sorted_bams_small_refset/*.bam
do
    echo $file
    sample_name=$(basename $file)
    sample_name=${sample_name%.*}
    samtools flagstat -O tsv $file > /nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/3_flagstat_small_refset/${sample_name}.tsv
done

for file in /nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/1_mapped_sorted_bams/*.bam
do
    echo $file
    sample_name=$(basename $file)
    sample_name=${sample_name%.*}
    samtools flagstat -O tsv $file > /nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2/3_flagstat/${sample_name}.tsv
done