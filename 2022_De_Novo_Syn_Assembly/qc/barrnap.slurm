#!/usr/bin/env bash
# 
#SBATCH --job-name=barrnap
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
#SBATCH --comment="Barrnap on de novo syn isolates"
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
#SBATCH -o logs/barrnap/%j_slurm_output.txt
#SBATCH -e logs/barrnap/%j_slurm_error.txt


wd="/nfs/chisholmlab001/kve/2022_De_Novo_Syn_Assembly/metagenomic_binning"

source activate barrnap

for bindir in $wd/*/*_bins/*;
do 
    binname=$(basename "$bindir")
    contigs=${bindir}/${binname}.fa
    rrnafasta=${bindir}/${binname}_rrna.fa
    rrnagff=${bindir}/${binname}_rrna.gff

    echo $contigs

    barrnap -o ${rrnafasta} < ${contigs} > ${rrnagff}

done
