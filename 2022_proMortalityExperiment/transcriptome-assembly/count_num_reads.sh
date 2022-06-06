

read_dir="/nfs/chisholmlab001/chisholmlab/experiment_repository/2019/ProMo_transcriptomics/210921_SONYA_SHEEAN_21_MARINEALGALCOMMUNITY_RNA_STRDPOLYA_40M_PE100_NOVASEQ-20220216T172330Z"

output="/home/kve/scripts/mini_projects/2022_proMortalityExperiment/transcriptome-assembly/counts.tsv"

for path in $read_dir/*
do
    filename=$(basename $path)
    filestem="${filename%%.*}"
    num_reads=$(zcat $path | echo $((`wc -l`/4)))
    echo -e "${filestem}\t${num_reads}" >> $output
done
