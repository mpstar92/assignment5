# scripts/variability.py (optimierte Version)
import pandas as pd
import numpy as np
from collections import Counter
from Bio import AlignIO

window_size = int(snakemake.config["window_size"])
alignment = AlignIO.read(snakemake.input[0], "fasta")

def shannon_entropy(column):
    counts = Counter(column)
    total = sum(counts.values())
    freqs = [count / total for count in counts.values()]
    return -sum(f * np.log2(f) for f in freqs if f > 0)

columns = list(zip(*[list(rec.seq) for rec in alignment]))

entropies = [shannon_entropy(col) for col in columns]

entropy_df = pd.DataFrame({
    "pos": range(len(entropies)),
    "entropy": entropies
})
entropy_df["avg_entropy"] = entropy_df["entropy"].rolling(window=window_size, min_periods=1).mean()

result = entropy_df.iloc[::window_size][["pos", "avg_entropy"]].rename(columns={"pos": "start", "avg_entropy": "entropy"})

result.to_csv(snakemake.output[0], sep="\t", index=False)