import pandas as pd
import numpy as np
from pathlib import Path
import gffpandas.gffpandas as gffpd
from Bio import SeqIO, pairwise2, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq3
from Bio.Seq import MutableSeq, Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from BCBio import GFF
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

# conda activate 2021_Sar11Pro_RNAseq_Project_Annotations

def df2gff(annotation_df, outpath, genome_record=None):
    gff_df = annotation_df[['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase']]
    attributes_df = annotation_df.drop(['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase'], axis=1)
    log.debug(f"gff columns: {gff_df.columns}")
    log.debug(f"attributes: columns {attributes_df.columns}")
    records = attributes_df.to_dict('records')
    gff_df['attributes'] = [';'.join(f"{k}={v}" for k,v in d.items()) for d in records]
    gff_df.to_csv(outpath, sep='\t', index=False, header=False)
    if genome_record:
        gff_df['seq_id'] = genome_record.id
        line_prepender(outpath, f"##sequence-region {genome_record.id} 1 {len(genome_record.seq)}")
    line_prepender(outpath, "##gff-version 3")

def gff2df(gff_path):
    annotation = gffpd.read_gff3(gff_path)
    attributes_df = annotation.attributes_to_columns()
    return attributes_df

def readFASTA(fasta_path):
    ref_genome_record = SeqIO.read(fasta_path, 'fasta')
    return ref_genome_record

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

def extract_cds_na(gff_df, ref_seq, fasta_outpath):
    assert isinstance(gff_df, pd.DataFrame)
    assert isinstance(ref_seq, Seq)
    gff_cols = gff_df.columns.to_list()
    assert "start" in gff_cols
    assert "end" in gff_cols
    assert "type" in gff_cols
    assert "CDS" in gff_df["type"].to_list()

    cds_df = gff_df[gff_df["type"] == 'CDS']
    gene_records = []

    for i, gene in cds_df.iterrows():
        na_seq = ref_seq[int(gene['start'])-1:int(gene['end'])]
        gene_seq_record = SeqRecord(na_seq, id=gene['ID'], name=str(gene['Name']), 
                                    description=f"{gene['start']}|{gene['strand']}|{gene['end']}")
        gene_records.append(gene_seq_record)

    SeqIO.write(gene_records, fasta_outpath, "fasta")


def extract_cds_aa(gff_df, ref_seq, fasta_outpath):
    assert isinstance(gff_df, pd.DataFrame)
    assert isinstance(ref_seq, Seq)
    gff_cols = gff_df.columns.to_list()
    assert "start" in gff_cols
    assert "end" in gff_cols
    assert "type" in gff_cols
    assert "CDS" in gff_df["type"].to_list()

    cds_df = gff_df[gff_df["type"] == 'CDS']
    gene_records = []

    for i, gene in cds_df.iterrows():
        na_seq = ref_seq[int(gene['start'])-1:int(gene['end'])]
        if gene['strand'] == '-':
            na_seq = na_seq.reverse_complement()
        aa_seq = na_seq.translate(stop_symbol="")
        gene_seq_record = SeqRecord(aa_seq, id=gene['ID'], name=str(gene['Name']), 
                                    description=f"{gene['start']}|{gene['strand']}|{gene['end']}")
        gene_records.append(gene_seq_record)

    SeqIO.write(gene_records, fasta_outpath, "fasta")
    

def run_blast_annotations(gff_df, ref_seq, database_outpath):
    assert isinstance(gff_df, pd.DataFrame)
    assert isinstance(ref_seq, Seq)
    gff_cols = gff_df.columns.to_list()
    assert "start" in gff_cols
    assert "end" in gff_cols
    assert "type" in gff_cols
    assert "CDS" in gff_df["type"].to_list()
    cds_df = gff_df[gff_df["type"] == 'CDS']
    data = []
    for i, gene in cds_df.iterrows():
        na_seq = ref_seq[int(gene['start'])-1:int(gene['end'])]
        if gene['strand'] == '-':
            na_seq = na_seq.reverse_complement()
        aa_seq = na_seq.translate(stop_symbol="")
        query_handle = NCBIWWW.qblast(program="blastp", database="nr", sequence=aa_seq) # TODO: consider "refseq_protein" database
        rec = NCBIXML.read(query_handle)
        gene_ID = gene['ID']
        best_reference = ""
        best_hit = ""
        best_e = np.nan
        best_score = np.nan
        for d in rec.descriptions:
            if 'hypothetical' not in d.title:
                best_reference = d.title.split("|")[1]
                best_hit = d.title.split("|")[2]
                best_e = d.e
                best_score = d.score
                break
        print(f"best ref: {best_reference}\t e: {best_e}\t score: {best_score}\t hit: {best_hit}\n\n")
        # header: ['gene_ID', 'best_hit', 'best_reference', 'best_e', 'best_score']
        data.append([gene_ID, best_hit, best_reference, best_e, best_score])
    df = pd.DataFrame(data, columns=['gene_ID', 'best_hit', 'best_reference', 'best_e', 'best_score'])
    df.to_csv(database_outpath, sep="\t")

def parse_blast_annotations(blast_xml, extracted_cds_fasta, table_outpath):

    with open(extracted_cds_fasta) as in_seq_handle:
        cds_records = list(SeqIO.parse(in_seq_handle, "fasta"))

    # with open(blast_xml, "r") as query_handle:
    #     blast_records = NCBIXML.parse(query_handle)
    #     len_blast_records = sum(1 for _ in blast_records)
    #     print(len_blast_records, len(cds_records))
    #     assert len_blast_records == len(cds_records)

    data = []
    with open(blast_xml, "r") as query_handle:
        blast_records = NCBIXML.parse(query_handle)
        for balst_rec, seq_rec in zip(blast_records, cds_records):

            gene_ID = seq_rec.id
            best_hit = ""
            best_hit_annot = ""
            best_e = np.nan
            best_score = np.nan
            best_pct_ident = np.nan

            for d in balst_rec.alignments:
                if 'hypothetical' not in d.hit_def:
                    best_hit = d.hit_id
                    best_hit_annot = d.hit_def
                    best_e = d.hsps[0].expect
                    best_score = d.hsps[0].score
                    best_pct_ident = d.hsps[0].identities / d.hsps[0].align_length
                    break
         
            row = [gene_ID, best_hit, best_hit_annot, best_e, best_score, best_pct_ident]
            print("\t".join([str(x) for x in row]))
            data.append(row)

    header = ['gene_ID', 'best_hit', 'best_hit_annot', 'best_e', 'best_score', 'best_pct_ident']
    df = pd.DataFrame(data, columns=header)
    df.to_csv(table_outpath, sep="\t", index=False)


def get_seq_record(fasta_path, gff_path):
    
    limit_info = dict(gff_type=["CDS"])

    in_seq_handle = open(fasta_path)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    record = None
    in_handle = open(gff_path)
    for rec in GFF.parse(in_handle, base_dict=seq_dict, limit_info=limit_info):
        record = rec
    in_handle.close()

    return record

def extract_Sar11Pro_project_CDS():
    genome_input_dir = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/input_data")

    genome_output_dir = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations")
    ncbi_MIT9301_resource_path = genome_output_dir / "MIT9301_Genomic_Resources" / "NCBI"
    ncbi_HTCC7211_resource_path = genome_output_dir / "HTCC7211_Genomic_Resources" / "NCBI"
    
    gff_path = ncbi_MIT9301_resource_path / "NCBI_MIT9301.gff"
    fasta_path = ncbi_MIT9301_resource_path / "NCBI_MIT9301.fna"
    fna_outpath = ncbi_MIT9301_resource_path / "MIT9301_extracted_cds.fna"
    faa_outpath = ncbi_MIT9301_resource_path / "MIT9301_extracted_cds.faa"

    gff_df = gff2df(gff_path)
    ref_seq = readFASTA(fasta_path)
    extract_cds_na(gff_df, ref_seq.seq, fna_outpath)
    extract_cds_aa(gff_df, ref_seq.seq, faa_outpath)

    gff_path = ncbi_HTCC7211_resource_path / "HTCC7211.gff"
    fasta_path = ncbi_HTCC7211_resource_path / "HTCC7211.fna"
    fna_outpath = ncbi_HTCC7211_resource_path / "HTCC7211_extracted_cds.fna"
    faa_outpath = ncbi_HTCC7211_resource_path / "HTCC7211_extracted_cds.faa"
    
    gff_df = gff2df(gff_path)
    ref_seq = readFASTA(fasta_path)
    extract_cds_na(gff_df, ref_seq.seq, fna_outpath)
    extract_cds_aa(gff_df, ref_seq.seq, faa_outpath)

def compare_gffs():
    improving_annotations_dir = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations")
    cluster_ref_path = improving_annotations_dir / "MIT9301_Genomic_Resources" / "Cluster_MIT9301_ref.fna"
    cluster_gff_path = improving_annotations_dir / "MIT9301_Genomic_Resources" / "Cluster_MIT9301.gff"
    ncbi_ref_path = improving_annotations_dir / "MIT9301_Genomic_Resources" / "NCBI_MIT9301_ref.fna"
    ncbi_gff_path = improving_annotations_dir / "MIT9301_Genomic_Resources" / "NCBI_MIT9301.gff"

    cluster_record = get_seq_record(cluster_ref_path, cluster_gff_path)
    ncbi_record = get_seq_record(ncbi_ref_path, ncbi_gff_path)

    cluster_hypotheticals = 0
    cluster_cds = 0
    for sf in cluster_record.features:
        cluster_cds += 1
        if 'hypothetical' in sf.qualifiers['product'][0]:
            cluster_hypotheticals += 1

    print(f'cluster hypotheticals: \t {cluster_hypotheticals}')
    print(f'cluster cds: \t\t {cluster_cds}')
    

    ncbi_hypotheticals = 0
    ncbi_cds = 0
    for i, sf in enumerate(ncbi_record.features):
        ncbi_cds += 1
        if sf.type == 'CDS':
            if 'hypothetical' in sf.qualifiers['product'][0]:
                ncbi_hypotheticals += 1
            
    print(f'ncbi hypotheticals: \t {ncbi_hypotheticals}')
    print(f'ncbi cds: \t\t {ncbi_cds}')

    cluster_annotation_df = gff2df(cluster_gff_path)
    ncbi_annotation_df = gff2df(ncbi_gff_path)

    cluster_annotation_df = cluster_annotation_df[cluster_annotation_df['type'] == 'CDS']
    ncbi_annotation_df = ncbi_annotation_df[ncbi_annotation_df['type'] == 'CDS']

    cluster_annotation_df.columns = pd.MultiIndex.from_product([['cluster'], cluster_annotation_df.columns])
    ncbi_annotation_df.columns = pd.MultiIndex.from_product([['ncbi'], ncbi_annotation_df.columns])

    cluster_annotation_df['shared', 'start_stop_ident'] = cluster_annotation_df['cluster']['start'].astype(str) + "|" + cluster_annotation_df['cluster']['strand'].astype(str) + "|" + cluster_annotation_df['cluster']['end'].astype(str)
    ncbi_annotation_df['shared', 'start_stop_ident'] = ncbi_annotation_df['ncbi']['start'].astype(str) + "|" + ncbi_annotation_df['ncbi']['strand'].astype(str) + "|" + ncbi_annotation_df['ncbi']['end'].astype(str)

    joint_cds_df = pd.merge(cluster_annotation_df, ncbi_annotation_df, on=[('shared', 'start_stop_ident')], how='outer')
    joint_cds_df = joint_cds_df.reset_index()

    joint_cds_df['metadata', 'cluster_hypothetical'] = joint_cds_df['cluster']['product'].str.contains('hypothetical')
    joint_cds_df['metadata', 'ncbi_hypothetical'] = joint_cds_df['ncbi']['product'].str.contains('hypothetical')
    joint_cds_df['metadata', 'annotated_in_neither'] = joint_cds_df['metadata']['cluster_hypothetical'] & joint_cds_df['metadata']['ncbi_hypothetical']
    joint_cds_df['metadata', 'annotated_in_both'] = ~(joint_cds_df['metadata']['cluster_hypothetical'] | joint_cds_df['metadata']['ncbi_hypothetical'])

    joint_cds_df['metadata', 'in_cluster_gff'] = joint_cds_df['cluster']['product'].notna()
    joint_cds_df['metadata', 'in_ncbi_gff'] = joint_cds_df['ncbi']['product'].notna()
    joint_cds_df['metadata', 'in_both'] = joint_cds_df['metadata']['in_cluster_gff'] & joint_cds_df['metadata']['in_ncbi_gff']

    joint_cds_df.to_csv(improving_annotations_dir / "gff_comparison.tsv", sep='\t', index=False)

    joint_cds_df[joint_cds_df['metadata', 'annotated_in_neither'] == True].to_csv(improving_annotations_dir / "gff_comparison_all_hypothetical.tsv", sep='\t', index=False)
    joint_cds_df[joint_cds_df['metadata', 'annotated_in_both'] == True].xs(key='product', axis=1, level=1).to_csv(improving_annotations_dir / "gff_comparison_products.tsv", sep='\t', index=False)
    
    cluster_df_w_ncbi = joint_cds_df[joint_cds_df['cluster']['product'].notna()]
    cluster_df_w_ncbi = cluster_df_w_ncbi.set_index(("cluster", "ID"))

    cluster_df_w_ncbi.xs(key='product', axis=1, level=1).to_csv(improving_annotations_dir / "add_ncbi_gff_annotations_cluster_gff_index.tsv", sep='\t')

def save_blast_annotations():
    improving_annotations_dir = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations")
    ref_path = improving_annotations_dir / "MIT9301_Genomic_Resources" / "Cluster_MIT9301_ref.fna"
    gff_path = improving_annotations_dir / "MIT9301_Genomic_Resources" / "Cluster_MIT9301.gff"

    database_outpath = improving_annotations_dir / "python_blastp_nr_results.tsv"

    # qresults = SearchIO.parse(blast_results_path, 'blast-xml')
    # search_dict = SearchIO.to_dict(qresults)

    gff_df = gff2df(gff_path)
    ref_rec = readFASTA(ref_path)

    run_blast_annotations(gff_df, ref_rec.seq, database_outpath)

def read_blast_annotations():
    improving_annotations_dir = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations")
    extracted_cds_path = improving_annotations_dir / "MIT9301_extracted_cds.fna"

    blast_xml = improving_annotations_dir / "MIT9301_cds_blastp_uniprot.xml"
    database_outpath = improving_annotations_dir / "blastp_swissprot_results.tsv"
    parse_blast_annotations(blast_xml, extracted_cds_path, database_outpath)

    # this XML was corrupted and uncompleted so I am writing over it. 
    # blast_xml = improving_annotations_dir / "MIT9301_cds_blastn.xml"
    # database_outpath = improving_annotations_dir / "blastn_nt_results.tsv"
    # parse_blast_annotations(blast_xml, extracted_cds_path, database_outpath)

def annotate_joined_df():
    joined_df_path = Path("/home/kve/scripts/mini_projects/2021_Sar11Pro_RNAseq_Project/analyse_media/joined_df.tsv")
    improving_annotations_dir = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations")
    ncbi_gff_annotations = improving_annotations_dir / "add_ncbi_gff_annotations_cluster_gff_index.tsv"
    swissprot_blast_results = improving_annotations_dir / "blastp_swissprot_results.tsv"

    output_df_path = improving_annotations_dir / "joined_df_with_ncbi_annotations.tsv"

    joined_df = pd.read_csv(joined_df_path, sep="\t", index_col=0)
    ncbi_gff_df = pd.read_csv(ncbi_gff_annotations, sep="\t", index_col=0)
    swissprot_blast_df = pd.read_csv(swissprot_blast_results, sep="\t", index_col="gene_ID")
    swissprot_blast_df = swissprot_blast_df.add_prefix('swissprot_blast_')

    joined_df = joined_df.join(ncbi_gff_df)
    joined_df = joined_df.join(swissprot_blast_df)

    joined_df.to_csv(output_df_path, sep="\t")


def addAnnotTable2GFF(gff_path, annot_df, gff_outpath):
    gff_df = gff2df(gff_path)

    log.debug(f"gff df:\n{gff_df.iloc[2]}\n")
    log.debug(f"added annotations df:\n{annot_df.iloc[0]}\n")
    
    annotation_df = gff_df.join(annot_df, on='ID', rsuffix='_annot')
    
    log.debug(f"joined df:\n{gff_df.iloc[2]}\n")

    df2gff(annotation_df, gff_outpath)

def addeggNOGannot2NCBIgff():
    improving_annotations_dir = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/improving_annotations")
    mit9301_ncbi_dir = improving_annotations_dir / "MIT9301_Genomic_Resources" / "NCBI"
    mit9301_eggnog_dir = improving_annotations_dir / "MIT9301_Genomic_Resources" / "MIT9301_NCBI_proteins_eggNOG_orthologs"
    htcc7211_ncbi_dir = improving_annotations_dir / "HTCC7211_Genomic_Resources" / "NCBI"
    htcc7211_eggnog_dir = improving_annotations_dir / "HTCC7211_Genomic_Resources" / "HTCC7211_NCBI_proteins_eggNOG_orthologs"

    mit9301_ncbi_gff_path = mit9301_ncbi_dir / "NCBI_MIT9301.gff"
    mit9301_eggnog_table_path = mit9301_eggnog_dir / "MIT9301_NCBI_proteins_eggNOG_orthologs_annotations_table.tsv"
    mit9301_gff_outpath = mit9301_eggnog_dir / "MIT9301_NCBI_proteins_eggNOG_orthologs.gff"

    htcc7211_ncbi_gff_path = htcc7211_ncbi_dir / "HTCC7211.gff"
    htcc7211_eggnog_table_path = htcc7211_eggnog_dir / "HTCC7211_NCBI_proteins_eggNOG_orthologs_annotations_table.tsv"
    htcc7211_gff_outpath = htcc7211_eggnog_dir / "HTCC7211_NCBI_proteins_eggNOG_orthologs.gff"

    mit9301_annot_df = pd.read_csv(mit9301_eggnog_table_path, sep='\t', index_col='query')
    htcc7211_annot_df = pd.read_csv(htcc7211_eggnog_table_path, sep='\t', index_col='query')

    addAnnotTable2GFF(mit9301_ncbi_gff_path, mit9301_annot_df, mit9301_gff_outpath)
    addAnnotTable2GFF(htcc7211_ncbi_gff_path, htcc7211_annot_df, htcc7211_gff_outpath)

if __name__ == "__main__":
    # read_blast_annotations()
    # compare_gffs()
    # read_blast_annotations()
    # annotate_joined_df()
    addeggNOGannot2NCBIgff()