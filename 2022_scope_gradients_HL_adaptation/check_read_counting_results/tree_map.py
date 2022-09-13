import pandas as pd
import plotly.graph_objects as go
from pathlib import Path

class MetagenomicAnalysis:
    def __init__(self, read_count_path, flagstat_dir, htseq_dir):
        self.read_counts = pd.read_table(read_count_path, names=['sample', 'type', 'count'], index_col=['sample', 'type'])
        self.flagstat_dict = {f.stem : pd.read_table(f, names=['qc_passed_reads', 'qc_failed_reads', 'label'], index_col='label') for f in Path(flagstat_dir).iterdir() if f.suffix == '.tsv'}
        
        self.htseq_metadata_dict = {}
        self.htseq_count_dict = {}
        for f in Path(htseq_dir).iterdir():
            if f.suffix == '.tsv':
                sample = f.stem
                htseq_df = pd.read_table(f, header=0, names=['gene_id', 'paired_counts', 'unpaired_counts'])

                htseq_counts_df = htseq_df[~htseq_df['gene_id'].str.contains("__")]
                htseq_counts_df = htseq_counts_df.set_index('gene_id')

                htseq_metadata_df = htseq_df[htseq_df['gene_id'].str.contains("__")]
                htseq_metadata_df['gene_id'] = htseq_metadata_df['gene_id'].apply(lambda x: x[2:])
                htseq_metadata_df = htseq_metadata_df.set_index('gene_id')
                htseq_metadata_df.loc['read_counts'] = htseq_counts_df.sum(axis=0)

                self.htseq_count_dict[sample] = htseq_counts_df
                self.htseq_metadata_dict[sample] = htseq_metadata_df
        

    def makeFlagstatIcicleChart(self, sample):
        flagstat_df = self.flagstat_dict[sample]
        tree = [
            {
                'label' : 'Reads from fastq',
                'parent' : '',
                'value' : self.read_counts.loc[sample].sum().loc['count'],
            },
            {
                'label' : 'Aligned to refs',
                'parent' : 'Reads from fastq',
                'value' : int(flagstat_df.loc["total (QC-passed reads + QC-failed reads)", 'qc_passed_reads']),
            },
            {
                'label' : 'primary',
                'parent' : 'Aligned to refs',
                'value' : int(flagstat_df.loc["primary", 'qc_passed_reads']),
            },
            {
                'label' : 'secondary',
                'parent' : 'Aligned to refs',
                'value' : int(flagstat_df.loc["secondary", 'qc_passed_reads']),
            },
            {
                'label' : 'supplementary',
                'parent' : 'Aligned to refs',
                'value' : int(flagstat_df.loc["supplementary", 'qc_passed_reads']),
            },
            {
                'label' : 'unpaired in sequencing',
                'parent' : 'primary',
                'value' : int(flagstat_df.loc["primary", 'qc_passed_reads']) - int(flagstat_df.loc['paired in sequencing', 'qc_passed_reads']),
            },
            {
                'label' : 'paired in sequencing',
                'parent' : 'primary',
                'value' : int(flagstat_df.loc['paired in sequencing', 'qc_passed_reads']),
            },
            {
                'label' : 'with itself and mate mapped',
                'parent' : 'paired in sequencing',
                'value' : int(flagstat_df.loc['with itself and mate mapped', 'qc_passed_reads']),
            },
            {
                'label' : 'singletons',
                'parent' : 'paired in sequencing',
                'value' : int(flagstat_df.loc['singletons', 'qc_passed_reads']),
            },
            {
                'label' : 'properly paired and mapped',
                'parent' : 'with itself and mate mapped',
                'value' : int(flagstat_df.loc['properly paired', 'qc_passed_reads']),
            }
        ]
        print(tree)
        fig = go.Figure(go.Icicle(
            labels = [n['label'] for n in tree],
            values = [n['value'] for n in tree],
            parents = [n['parent'] for n in tree],
            branchvalues="total",
            root_color="lightgrey"
        ))
        return fig

    def makeHTseqMetadataIcicleChart(self, sample):
        flagstat_df = self.flagstat_dict[sample]
        htseq_metadata_df = self.htseq_metadata_dict[sample]
        tree = [
            {
                'label' : 'Reads from fastq',
                'parent' : '',
                'value' : self.read_counts.loc[sample].sum().loc['count'],
            },
            {
                'label' : 'Aligned to refs',
                'parent' : 'Reads from fastq',
                'value' : int(flagstat_df.loc["total (QC-passed reads + QC-failed reads)", 'qc_passed_reads']),
            },
            {
                'label' : 'Paired',
                'parent' : 'Aligned to refs',
                'value' : int(htseq_metadata_df['paired_counts'].sum()),
            },
            {
                'label' : 'Unpaired',
                'parent' : 'Aligned to refs',
                'value' : int(htseq_metadata_df['unpaired_counts'].sum()),
            },
        ]

        for datum in htseq_metadata_df.index:
            tree.append({
                'label' : datum,
                'parent' : 'Paired',
                'value' : int(htseq_metadata_df.loc[datum, 'paired_counts']),
            })
            tree.append({
                'label' : datum,
                'parent' : 'Unpaired',
                'value' : int(htseq_metadata_df.loc[datum, 'unpaired_counts']),
            })
        print(tree)
        fig = go.Figure(go.Icicle(
            labels = [n['label'] for n in tree],
            values = [n['value'] for n in tree],
            parents = [n['parent'] for n in tree],
            branchvalues="total",
            root_color="lightgrey"
        ))
        return fig


wd = Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/intermediates_bowtie2")

analysis = MetagenomicAnalysis(wd / "3_flagstat" / "read_counts.tsv", wd / "3_flagstat" / "flagstat" / "all", wd / '2_htseq_concat_output_0_mapq')
analysis.makeHTseqMetadataIcicleChart('G1St3100m_S30').write_html("HTseqMetadataIcicleChart_G1St3100m_S30_MAPQ0.html")

analysis = MetagenomicAnalysis(wd / "3_flagstat" / "read_counts.tsv", wd / "3_flagstat" / "flagstat" / "all", wd / '2_htseq_concat_output_10_mapq')
analysis.makeHTseqMetadataIcicleChart('G1St3100m_S30').write_html("HTseqMetadataIcicleChart_G1St3100m_S30_MAPQ10.html")

analysis = MetagenomicAnalysis(wd / "3_flagstat" / "read_counts.tsv", wd / "3_flagstat" / "flagstat" / "all", wd / '2_htseq_concat_output_30_mapq')
analysis.makeHTseqMetadataIcicleChart('G1St3100m_S30').write_html("HTseqMetadataIcicleChart_G1St3100m_S30_MAPQ30.html")