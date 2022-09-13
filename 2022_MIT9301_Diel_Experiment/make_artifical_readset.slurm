#!/usr/bin/env bash
# 
#SBATCH --job-name=art
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 4-0
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#SBATCH -n 1
#
#
#SBATCH -o logs/art/%j_slurm_output.txt
#SBATCH -e logs/art/%j_slurm_error.txt

internal_standards_fasta="/nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/internal_standards.fa"
MIT9301_fasta="/nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/MIT9301_seq_db/MIT9301_IMG_2623620961.fna"

OUT_DIR="/nfs/chisholmlab001/kve/2022_MIT9301_Diel_Experiment/internal_standard_check/artifical_read_mapping"

# /home/kve/bin/art_bin_MountRainier/art_illumina -i ${internal_standards_fasta} -ss "HS20" -c 50000000 -l 40 --paired --mflen 325 -s 30 -o ${OUT_DIR}/"hiseq2000"/"art_internal_standards_"
/home/kve/bin/art_bin_MountRainier/art_illumina -i ${internal_standards_fasta} -ss "NS50" -c 50000000 -l 75 --paired --mflen 325 -s 30 -o ${OUT_DIR}/"nextseq500"/"art_internal_standards_"
# /home/kve/bin/art_bin_MountRainier/art_illumina -i ${MIT9301_fasta} -ss "HS20" -c 50000000 -l 40 --paired --mflen 325 -s 30 -o ${OUT_DIR}/"hiseq2000"/"art_MIT9301_"
/home/kve/bin/art_bin_MountRainier/art_illumina -i ${MIT9301_fasta} -ss "NS50" -c 50000000 -l 75 --paired --mflen 325 -s 30 -o ${OUT_DIR}/"nextseq500"/"art_MIT9301_"

