import pandas as pd
import numpy as np
from pathlib import Path
import gffpandas.gffpandas as gffpd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

wd = Path("/nfs/chisholmlab001/kve/2022_WH7803_variant_calling/genomics")
gff_path =  wd / "pacbio_WH7803x_sorted.gff"
ref_path = wd / "pacbio_WH7803x_sorted.fasta"
cds_aa_path = wd / "pacbio_WH7803x_sorted.faa" 



def create_attr_dict(attr_row):
    logging.debug(attr_row)
    d = {}
    for key_value_pair in attr_row.split(";"):
        k_v_list = key_value_pair.split(sep="=", maxsplit=1)
        if len(k_v_list) == 2:
            k, v = k_v_list
            d[k] = v
    return d

def generate_gff_df(gff_file):
    ## This method is heavily influenced from GFFpandas
    df = pd.read_csv(gff_file, sep="\t", comment="#",
            names=[
                "seq_id",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
    )
    logging.debug(f"gff_file: {gff_file}")
    logging.debug(f"df.columns: {df.columns}")
    logging.debug(f"df.attributes: {df.attributes}")
    attr_dict_series = df.attributes.apply(
        lambda attributes: create_attr_dict(attributes)
    )
    key_set = set()
    attr_dict_series.apply(
        lambda at_dic: key_set.update(list(at_dic.keys()))
    )
    for attr in sorted(list(key_set)):
        df[attr] = attr_dict_series.apply(
            lambda attr_dict: attr_dict.get(attr)
        )
        
    return df
    

annotation = generate_gff_df(gff_path).set_index(["seq_id", "ID"])

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