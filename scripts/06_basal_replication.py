"""Step 06 - Independent replication: the same contrast in basal cells.

Runs the identical lactating vs pregnant comparison on the basal cell
population of the same study. Two independent cell types showing the same
lactation signature is much stronger evidence than either 2-replicate
contrast alone.

Outputs:
  results/tables/de_results_basal_lactating_vs_pregnant_full.csv
  results/tables/de_genes_basal_significant.csv
"""
import gzip
from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv(gzip.open(ROOT / "data" / "GSE60450_Lactation-GenewiseCounts.txt.gz", "rt"), sep="\t")
counts = counts.rename(columns={c: c.split("_")[0] for c in counts.columns if c.startswith("MCL1")})
meta = pd.read_csv(ROOT / "data" / "processed" / "sample_metadata.csv")

basal = meta[meta["cell_type"] == "basal"]["sample"].tolist()
basal_counts = counts.set_index("EntrezGeneID")[basal]
keep = (basal_counts >= 10).sum(axis=1) >= 2
basal_counts = basal_counts[keep]
print(f"basal genes after filtering: {len(basal_counts):,}")

bmeta = meta[meta["cell_type"] == "basal"].set_index("sample").loc[basal]
dds = DeseqDataSet(counts=basal_counts.T, metadata=bmeta[["stage"]], design="~stage", refit_cooks=True)
dds.deseq2()
ds = DeseqStats(dds, contrast=["stage", "lactating", "pregnant"])
ds.summary()
res = ds.results_df.sort_values("padj")

gi = pd.read_csv(gzip.open(ROOT / "data" / "Mus_musculus.gene_info.gz", "rt"),
                 sep="\t", usecols=["GeneID", "Symbol"], dtype={"GeneID": str})
sym = gi.drop_duplicates("GeneID").set_index("GeneID")["Symbol"]
res.insert(0, "symbol", [sym.get(str(g), "") for g in res.index])

res.to_csv(TABLES / "de_results_basal_lactating_vs_pregnant_full.csv")
sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)]
sig.to_csv(TABLES / "de_genes_basal_significant.csv")
print(f"basal significant (padj<0.05, |log2FC|>=1): {len(sig):,} "
      f"(up: {(sig['log2FoldChange'] > 0).sum():,}, down: {(sig['log2FoldChange'] < 0).sum():,})")
