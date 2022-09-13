#!/usr/bin/env bash
# 
#SBATCH --job-name=blastn
#
# Number of tasks/cores for job
#SBATCH -n 20
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 24:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH --mem=250G
#
#SBATCH -o logs/blastn/%j_slurm_output.txt
#SBATCH -e logs/blastn/%j_slurm_error.txt

source activate blast

# cd /nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/het_seq_db
# database="all_het_seqs"
# output="internal_standards_het_blastn.tsv"

# blastn -query "/nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/internal_standards.fa" -db ${database} -evalue 10 -word_size 8 -outfmt 6 -out ${output}

# cd /nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/pro_syn_seq_db
# database="all_pro_syn_seqs"
# output="internal_standards_pro_syn_blastn.tsv"

# blastn -query "/nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/internal_standards.fa" -db ${database} -evalue 10 -word_size 8 -outfmt 6 -out ${output}


# cd /nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/NB0046_comparison_db
# database="nb0046"
# output="internal_standards_nb0046_blastn.tsv"

# blastn -query "/nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/internal_standards.fa" -db ${database} -evalue 10 -word_size 8 -outfmt 6 -out ${output}


# cd /nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/NB0046_comparison_genes_db
# database="nb0046_genes"
# output="internal_standards_nb0046_genes_blastn.tsv"

# blastn -query "/nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/internal_standards.fa" -db ${database} -evalue 10 -word_size 8 -outfmt 6 -out ${output}

cd /nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/MIT9301_seq_db
database="MIT9301_IMG_2623620961_seq"
output="internal_standards_MIT9301_seq_blastn.tsv"

blastn -query "/nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/internal_standards.fa" -db ${database} -evalue 10 -word_size 8 -outfmt 6 -out ${output}
