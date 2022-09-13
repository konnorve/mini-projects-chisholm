import pandas as pd
from pathlib import Path

seq_ids = []
dfs = []

for gff in Path("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/chosen_reference_gffs/").glob("*.gff"):
    df = pd.read_table(gff, comment="#", header=None, names=[
                "seq_id",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ])
    seq_ids.append(set(df.seq_id.to_list()))
    dfs.append(df)

a = len(set.union(*seq_ids))
b = sum([len(s) for s in seq_ids])
assert a == b
print(a, b)

file_contents = pd.concat(dfs, ignore_index=True).to_csv(sep="\t", header=False, index=False)

with open("/nfs/chisholmlab001/kve/2022_scope_gradients_HL_adaptation/inputs/reference_database/concat_reference_annotations.gff", 'w') as f:
    f.write("##gff-version 3\n"+file_contents)