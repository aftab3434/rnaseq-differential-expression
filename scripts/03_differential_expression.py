"""Step 03 - Differential expression with PyDESeq2.

Compares lactating vs pregnant luminal cells (2 vs 2 replicates).
PyDESeq2 is a Python re-implementation of DESeq2 (Muzellec et al. 2023,
Bioinformatics), matching DESeq2 v1.34 default behaviour. Small numeric
differences from R DESeq2 are expected and documented by its authors.

Outputs:
  results/tables/de_results_lactating_vs_pregnant_full.csv  - all tested genes
  results/tables/de_genes_significant.csv                   - padj < 0.05 and |log2FC| >= 1
  data/processed/normed_counts_luminal.csv                  - normalized counts for plotting
"""
import gzip
from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parent.parent
(ROOT / "results" / "tables").mkdir(parents=True, exist_ok=True)

counts = pd.read_csv(ROOT / "data" / "processed" / "counts_luminal_filtered.csv", index_col=0)
meta = pd.read_csv(ROOT / "data" / "processed" / "sample_metadata.csv")
meta = meta[meta["cell_type"] == "luminal"].set_index("sample").loc[counts.columns]

# PyDESeq2 wants samples as rows, genes as columns
dds = DeseqDataSet(
    counts=counts.T,
    metadata=meta[["stage"]],
    design="~stage",
    refit_cooks=True,
)
dds.deseq2()

# normalized counts for downstream plots
norm = pd.DataFrame(dds.layers["normed_counts"], index=counts.columns, columns=counts.index)
norm.T.to_csv(ROOT / "data" / "processed" / "normed_counts_luminal.csv")

ds = DeseqStats(dds, contrast=["stage", "lactating", "pregnant"])
ds.summary()
res = ds.results_df.sort_values("padj")

# Entrez GeneID -> gene symbol (NCBI gene_info)
gi = pd.read_csv(
    gzip.open(ROOT / "data" / "Mus_musculus.gene_info.gz", "rt"),
    sep="\t", usecols=["GeneID", "Symbol"], dtype={"GeneID": str},
)
sym = gi.drop_duplicates("GeneID").set_index("GeneID")["Symbol"]
res.insert(0, "symbol", [sym.get(str(g), "") for g in res.index])

res.to_csv(ROOT / "results" / "tables" / "de_results_lactating_vs_pregnant_full.csv")

sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)]
sig.to_csv(ROOT / "results" / "tables" / "de_genes_significant.csv")

up = (sig["log2FoldChange"] > 0).sum()
down = (sig["log2FoldChange"] < 0).sum()
print(f"\ntested genes: {res['padj'].notna().sum():,}")
print(f"significant (padj<0.05, |log2FC|>=1): {len(sig):,}  (up in lactation: {up:,}, down: {down:,})")
print("\ntop 10 by adjusted p-value:")
print(sig.head(10)[["symbol", "log2FoldChange", "padj"]].to_string())
