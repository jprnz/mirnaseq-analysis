# logging
sys.stderr = open(snakemake.log[0], "w")

import sys
import pandas as pd

from pathlib import Path

def read_file(file):
    return pd.read_table(
        file,
        index_col=0,
        header=None,
        names=["mirna_id", Path(file).stem],
        dtype=object)


# Get counts, combine and write
cnt = [read_file(v) for v in snakemake.input]
dat = pd.concat(cnt, axis=1)
dat.index.name = "mirna_id"
dat.to_csv(snakemake.output[0])
