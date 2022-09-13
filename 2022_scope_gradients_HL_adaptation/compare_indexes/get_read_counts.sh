
rm read_counts.tsv
for file in /nobackup1/chisholmlab/kve/2022_scope_gradients_HL_adaptation/scratch/1_trimmed_reads/*.fastq.gz
do
    file_name=$(basename $file)
    echo $file_name
    file_stem=${file_name%.fastq.gz}
    IFS='_' read -r -a a <<< $file_stem
    read_count=$(echo $(zcat $file|wc -l)/4|bc)
    echo -e "${a[0]}_${a[1]}\t${a[2]}_${a[3]}\t${read_count}" >> read_counts.tsv
done