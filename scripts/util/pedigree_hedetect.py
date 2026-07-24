import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "paf", help="PAF file with mapping to X5A."
)
args = parser.parse_args()

cols = ["qname","qlen","qstart","qend","strand",
        "tname","tlen","tstart","tend","nmatch","alen","mapq"]

paf = pd.read_csv(
    args.paf, sep = "\t",
    usecols = range(12),
    names = cols
    )

paf["identity"] = paf["nmatch"] / paf["alen"]

window = 100_000
paf["window"] = (paf["qstart"] // window) * window

# For each window, find the best-scoring target chromosome.
best = (
    paf.sort_values("nmatch", ascending=False)
        .groupby("window")
        .agg(
            best_chrom=("tname","first"),
            best_identity=("identity","first")
        )
        .reset_index()
    )

best["is_HE_candidate"] = best["best_chrom"].str.contains("C0[0-9]|C10")