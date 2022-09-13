import pandas as pd
from Bio import SeqIO
import re

# # duplicated region:
# repeat_file = "/nfs/chisholmlab001/kve/2022_Population_Genomics_9313/genome_duplication_analysis/potential_duplication.repeats"
# fasta_file = "/nfs/chisholmlab001/kve/2022_Population_Genomics_9313/genome_duplication_analysis/potential_duplication.fasta"
# gff_outfile = "/nfs/chisholmlab001/kve/2022_Population_Genomics_9313/genome_duplication_analysis/potential_duplication.gff"

# whole genome:
repeat_file = "/nfs/chisholmlab001/kve/2022_Population_Genomics_9313/genome_duplication_analysis/original_genome.repeats"
fasta_file = "/nfs/chisholmlab001/kve/genomic_resources/strains/pro/MIT9313/MIT9313.fna"
gff_outfile = "/nfs/chisholmlab001/kve/2022_Population_Genomics_9313/genome_duplication_analysis/original_genome.gff"

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

def df2gff3(annotation_df, outpath):
    if 'attributes' not in annotation_df.columns:
        info_df = annotation_df.drop(['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase'], axis=1)
        annotation_df['attributes'] = info_df.apply(lambda r: ";".join([f"{k}={v}" for k, v in r.items()]), axis=1)
    gff_df = annotation_df[['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']]
    gff_df.to_csv(outpath, sep='\t', index=False, header=False)
    line_prepender(outpath, "##gff-version 3")
    

def parse_repeat_file(repeat_file):
    with open(repeat_file) as f:
        nl = next(f, None)
        while nl is not None:
            pl = nl
            nl = next(f, None)
            if pl.strip() == 'Exhaustive Exact Matches:':
                header = nl.split()
                nl = next(f, None)
                matches = []
                while nl is not None:
                    if 'r' in nl:
                        split_line = re.sub("r", "", nl).split()
                        if len(split_line) == 3:
                            match = dict(zip(['start1', 'end2', 'length'], map(int, split_line)))
                            match['direction'] = 'rev_complement'
                            match['end1'] = match['start1'] + match['length']
                            match['end2'] = match['end2'] + 1
                            match['start2'] = match['end2'] - match['length']
                            match['strand1'] = '+'
                            match['strand2'] = '-'
                            matches.append(match)
                    else:
                        split_line = nl.split()
                        if len(split_line) == 3:
                            match = dict(zip(['start1', 'start2', 'length'], map(int, split_line)))
                            match['direction'] = 'regular'
                            match['end1'] = match['start1'] + match['length']
                            match['end2'] = match['start2'] + match['length']
                            match['strand1'] = '+'
                            match['strand2'] = '+'
                            matches.append(match)
                    nl = next(f, None)
                return pd.DataFrame(matches)

df = parse_repeat_file(repeat_file)
seqrecord = list(SeqIO.parse(fasta_file, "fasta"))[0]

df['end1'] = df.start1+df.length

df['Seq1'] = df.apply(lambda r: seqrecord.seq[r['start1']-1:r['end1']-1], axis=1)
df['Seq2'] = df.apply(lambda r: seqrecord.seq[r['start2']-1:r['end2']-1], axis=1)
df['Seq2'] = df.apply(lambda r: r['Seq2'].reverse_complement() if r['direction']=='rev_complement' else r['Seq2'], axis=1)

df['Match'] = df.apply(lambda r: r['Seq1']==r['Seq2'], axis=1)

print(df.shape)

assert len(df[df.Match == False]) == 0

df['ID'] = [f"repeat_{i}" for i in range(len(df))]


df = df.query('280000 <= start1 <= 300000 | 280000 <= start2 <= 300000')


start1_df = pd.concat([df.start1, df.end1, df.strand1, df.ID, df.direction], axis=1)
start2_df = pd.concat([df.start2, df.end2, df.strand2, df.ID, df.direction], axis=1)

start1_df.columns = ['start', 'end', 'strand', 'ID', 'direction']
start2_df.columns = ['start', 'end', 'strand', 'ID', 'direction']

repeats_df = pd.concat([start1_df, start2_df])

repeats_df['seq_id'] = seqrecord.id
repeats_df['source'] = 'mummer_repeat'
repeats_df['type'] = 'repeat'
repeats_df['score'] = '.'
repeats_df['phase'] = '.'

print(repeats_df)

df2gff3(repeats_df, gff_outfile)