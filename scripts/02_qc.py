"""Step 02 - Quality control and low-count filtering.

- Reads the raw count matrix and attaches sample metadata
  (cell type and developmental stage, from the GEO series matrix order).
- Plots library sizes (total reads per sample).
- Plots the count distribution before and after filtering.
- Filters genes with very low counts (standard DESeq2 practice:
  keep genes with >= 10 counts in at least 2 samples).
- Writes the luminal-cell subset used for differential expression:
  data/processed/counts_luminal_filtered.csv and data/processed/sample_metadata.csv
"""
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA, FIGS = ROOT / "data", ROOT / "results" / "figures"
(ROOT / "data" / "processed").mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv(gzip.open(DATA / "GSE60450_Lactation-GenewiseCounts.txt.gz", "rt"), sep="\t")
counts = counts.rename(columns={c: c.split("_")[0] for c in counts.columns if c.startswith("MCL1")})

# Sample order in the file matches the GEO series matrix:
# DG..DL = luminal (virgin 1-2, pregnant 1-2, lactating 1-2)
# LA..LF = basal  (virgin 1-2, pregnant 1-2, lactating 1-2)
stages = ["virgin", "virgin", "pregnant", "pregnant", "lactating", "lactating"]
meta = pd.DataFrame({
    "sample": [f"MCL1-D{x}" for x in "GHIJKL"] + [f"MCL1-L{x}" for x in "ABCDEF"],
    "cell_type": ["luminal"] * 6 + ["basal"] * 6,
    "stage": stages * 2,
})
sample_cols = meta["sample"].tolist()
assert set(sample_cols) <= set(counts.columns), "sample columns not found in count file"

# --- library sizes ---
lib_sizes = counts[sample_cols].sum() / 1e6
fig, ax = plt.subplots(figsize=(8, 4))
colors = meta["cell_type"].map({"luminal": "#1f77b4", "basal": "#d62728"})
ax.bar(range(len(lib_sizes)), lib_sizes.values, color=colors)
ax.set_xticks(range(len(lib_sizes)))
ax.set_xticklabels([f"{s}\n{t}" for s, t in zip(meta["sample"], meta["stage"])], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Total counts (millions)")
ax.set_title("Library size per sample (blue = luminal, red = basal)")
fig.tight_layout()
fig.savefig(FIGS / "library_sizes.png", dpi=150)
plt.close(fig)
print(f"library sizes (millions):\n{lib_sizes.round(1).to_string()}")

# --- low-count filtering (on the luminal subset used for DE) ---
lum = meta[meta["cell_type"] == "luminal"]["sample"].tolist()
lum_counts = counts.set_index("EntrezGeneID")[lum]
keep = (lum_counts >= 10).sum(axis=1) >= 2
print(f"\ngenes before filtering: {len(lum_counts):,}")
print(f"genes after filtering (>=10 counts in >=2 samples): {keep.sum():,}")

fig, ax = plt.subplots(figsize=(7, 4))
for label, sub in [("before filtering", lum_counts), ("after filtering", lum_counts[keep])]:
    vals = np.log10(sub.values.flatten() + 1)
    ax.hist(vals, bins=80, histtype="step", label=label, linewidth=1.5)
ax.set_xlabel("log10(count + 1)")
ax.set_ylabel("frequency")
ax.set_title("Count distribution, luminal samples")
ax.legend()
fig.tight_layout()
fig.savefig(FIGS / "count_distribution.png", dpi=150)
plt.close(fig)

lum_counts[keep].to_csv(ROOT / "data" / "processed" / "counts_luminal_filtered.csv")
meta.to_csv(ROOT / "data" / "processed" / "sample_metadata.csv", index=False)
print("\nwrote data/processed/counts_luminal_filtered.csv and sample_metadata.csv")
