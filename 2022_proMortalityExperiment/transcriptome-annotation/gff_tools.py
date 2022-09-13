import pandas as pd
import numpy as np
from pathlib import Path
import gffpandas.gffpandas as gffpd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

wd = Path("/nfs/chisholmlab001/kve/2022_proMortalityExperiment/transcriptome_assembly/transcriptome_annotation")
gff_path =  wd / "eggNOGmapper" / "nr_transcriptome.fasta.transdecoder.CDS.gff3"
ref_path = wd / "finished_product" / "nr_transcriptome.fasta"
cds_aa_path = wd / "eggNOGmapper" / "nr_transcriptome_cds.faa" 

annotation = gffpd.read_gff3(gff_path).attributes_to_columns().set_index(["seq_id", "ID"])

aa_records = []
with open(ref_path) as handle:
    for transcript_record in SeqIO.parse(handle, "fasta"):
        if transcript_record.id in annotation.index.get_level_values(0):
            transcript_annotations = annotation.loc[transcript_record.id]
            for i in transcript_annotations.index.values:
                print(transcript_record.id, i, sep="\t")
                gene_annotations = transcript_annotations.loc[i]
                na_seq = transcript_record.seq[gene_annotations['start']-1:gene_annotations['end']]
                if gene_annotations['strand'] == '-':
                    na_seq = na_seq.reverse_complement()
                aa_seq = na_seq.translate(stop_symbol="")
                gene_record = SeqRecord(aa_seq, id=i)
                aa_records.append(gene_record)

with open(cds_aa_path, 'w') as handle:
    SeqIO.write(aa_records, handle, "fasta")