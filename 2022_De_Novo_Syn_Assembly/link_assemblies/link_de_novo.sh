
# link_dir=/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/all_syn_fastas
link_dir=/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/all_de_novo_bins

for d in /nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/metagenomic_binning/*
do 
    genome_names=$(grep  "Cyanobacteria" $d/gtdbtk_classify_wf/gtdbtk.bac120.summary.tsv | cut -f1)
    echo $d
    for genome_name in $genome_names
    do
        ln -s ${d}/bin_fastas/${genome_name}.fa ${link_dir}
    done
done