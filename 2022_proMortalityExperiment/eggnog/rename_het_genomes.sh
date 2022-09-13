
COUNTER=1
for f in /nfs/chisholmlab001/skearney/201111Alm/DAS_OUT/DAS_OUT_DASTool_bins/*.fa
do
    name=$(basename $f)
    name=${name%.*}
    new_name=$(printf "%02d" $COUNTER)
    new_name="bin_$new_name"

    echo $COUNTER
    echo $name
    echo $new_name

    cp $f "/nfs/chisholmlab001/kve/2022_proMortalityExperiment/input_data_ribozero/het_bins/$new_name.fna"
    cp "/nfs/chisholmlab001/skearney/201111Alm/DAS_OUT/DAS_OUT_DASTool_bins/checkm_out/bins/$name/genes.gff" "/nfs/chisholmlab001/kve/2022_proMortalityExperiment/input_data_ribozero/het_bins/$new_name.gff"
    cp "/nfs/chisholmlab001/skearney/201111Alm/DAS_OUT/DAS_OUT_DASTool_bins/checkm_out/bins/$name/genes.faa" "/nfs/chisholmlab001/kve/2022_proMortalityExperiment/input_data_ribozero/het_bins/$new_name.faa"

    COUNTER=$[$COUNTER +1]
done
