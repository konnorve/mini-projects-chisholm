import re
from pathlib import Path
import os

target_regex = re.compile(r'SS\d\d\d\d_(S\d+)_L002_(R1|R2)_001.fastq.gz')
target_regex = re.compile(r'2021_Sar11Pro_RNAseq_Project_(\d\d)_(1|2)_sequence.fastq.gz')

target_dir = Path('/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/input_data/raw_reads')

for f in target_dir.iterdir():
    so = target_regex.search(f.name)
    new_name = f"S{so.group(1)}_R{so.group(2)}.fastq.gz"
    os.rename(target_dir / f.name, target_dir / new_name)