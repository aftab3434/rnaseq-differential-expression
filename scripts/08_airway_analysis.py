#!/usr/bin/env python3
"""Script 08: Independent second-dataset analysis - the canonical `airway`
benchmark (Himes et al. 2014, GSE52778; counts via recount3, study SRP033351).

Design (matching the DESeq2 vignette): ~ cell_line + treatment, contrasting
dexamethasone (Dex) vs untreated in 4 human airway smooth muscle cell lines
(4 treated vs 4 untreated replicates).

Why this dataset: it is the standard DESeq2 teaching dataset and, unlike the
primary mouse dataset (n=2 per group), it has 4 replicates per condition,
giving DESeq2 a much better basis for dispersion estimation.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

RANDOM_SEED = 42
SYMBOLS = {"ENSG00000165995": "CACNB2", "ENSG00000146250": "PRSS35",
           "ENSG00000152583": "SPARCL1", "ENSG00000101347": "SAMHD1",
           "ENSG00000120129": "DUSP1", "ENSG00000189221": "MAOA",
           "ENSG00000157150": "TIMP4", "ENSG00000211445": "GPX3",
           "ENSG00000162692": "VCAM1", "ENSG00000196136": "SERPINA3"}
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 1.0
DATA = "data/airway_gene_sums.tsv.gz"
OUT = "results/airway"

# GEO series-matrix order: for each cell line: untreated, Dex, Alb, Alb_Dex
CELL_LINES = ["N61311", "N052611", "N080611", "N061011"]
TREATMENTS = ["untrt", "trt", "alb", "alb_trt"]  # trt = dexamethasone

print("Loading recount3 gene counts ...")
counts = pd.read_csv(DATA, sep="\t", comment="#", index_col=0)
counts = counts.T  # samples x genes
counts.index.name = "run"

meta = pd.DataFrame(
    [(run, CELL_LINES[i // 4], TREATMENTS[i % 4]) for i, run in enumerate(counts.index)],
    columns=["run", "cell_line", "treatment"],
).set_index("run")

# Canonical airway subset: untreated vs Dex only (8 samples, 4 vs 4)
keep = meta["treatment"].isin(["untrt", "trt"])
counts, meta = counts.loc[keep], meta.loc[keep]
print(f"Subset: {counts.shape[0]} samples x {counts.shape[1]} genes")
print(meta[["cell_line", "treatment"]])

# Filter low counts (DESeq2 vignette rule: row sum >= 10)
counts = counts.loc[:, counts.sum(axis=0) >= 10]
print(f"Genes after filtering: {counts.shape[1]}")

dds = DeseqDataSet(
    counts=counts.astype(int),
    metadata=meta,
    design="~ cell_line + treatment",
    ref_level=["treatment", "untrt"],
    refit_cooks=True,
)
dds.deseq2()

ds = DeseqStats(dds, contrast=["treatment", "trt", "untrt"])
ds.summary()
res = ds.results_df.sort_values("padj")
res.to_csv(f"{OUT}/airway_dex_vs_untreated_full.csv")

sig = res[(res.padj < PADJ_CUTOFF) & (res.log2FoldChange.abs() >= LFC_CUTOFF)]
sig_fdr_only = res[res.padj < PADJ_CUTOFF]
sig.to_csv(f"{OUT}/airway_significant_genes.csv")
print(f"\nFDR<{PADJ_CUTOFF}: {len(sig_fdr_only)} | + |L2FC|>=1: {len(sig)} "
      f"(up {int((sig.log2FoldChange>0).sum())}, down {int((sig.log2FoldChange<0).sum())})")
print("\nTop 10 by adjusted p-value:")
print(res.head(10)[["baseMean", "log2FoldChange", "pvalue", "padj"]])

# QC: library sizes
lib = counts.sum(axis=1)
fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(range(len(lib)), lib.values / 1e6,
       color=["#4C72B0" if t == "untrt" else "#C44E52" for t in meta["treatment"]])
ax.set_xticks(range(len(lib)))
ax.set_xticklabels([f"{cl}\n{t}" for cl, t in zip(meta.cell_line, meta.treatment)],
                   rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Total counts (millions)")
ax.set_title("airway: library sizes (blue=untreated, red=Dex)")
plt.tight_layout(); plt.savefig(f"{OUT}/library_sizes.png", dpi=150); plt.close()

# PCA on log2 normalized counts - expect cell line to dominate, treatment secondary
from sklearn.decomposition import PCA
norm = np.asarray(dds.layers["normed_counts"])
X = np.log2(norm + 1.0)
X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-8)
pca = PCA(n_components=2, random_state=RANDOM_SEED)
coords = pca.fit_transform(X)
fig, ax = plt.subplots(figsize=(6, 5))
for cl in CELL_LINES:
    m = meta.cell_line == cl
    for t, mk in [("untrt", "o"), ("trt", "s")]:
        mm = m & (meta.treatment == t)
        ax.scatter(coords[mm, 0], coords[mm, 1], marker=mk, s=90,
                   label=f"{cl} {t}")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("airway: PCA (circle=untreated, square=Dex)")
ax.legend(fontsize=7, ncol=2)
plt.tight_layout(); plt.savefig(f"{OUT}/pca_plot.png", dpi=150); plt.close()

# Volcano
fig, ax = plt.subplots(figsize=(7, 6))
plot_df = res.dropna(subset=["padj"])
colors = np.where(plot_df.padj < PADJ_CUTOFF,
         np.where(plot_df.log2FoldChange >= LFC_CUTOFF, "#C44E52",
         np.where(plot_df.log2FoldChange <= -LFC_CUTOFF, "#4C72B0", "grey")), "lightgrey")
ax.scatter(plot_df.log2FoldChange, -np.log10(plot_df.padj), c=colors, s=6, alpha=0.6)
ax.axhline(-np.log10(PADJ_CUTOFF), ls="--", lw=0.8, c="k")
ax.axvline(LFC_CUTOFF, ls="--", lw=0.8, c="k"); ax.axvline(-LFC_CUTOFF, ls="--", lw=0.8, c="k")
SYMBOLS = {"ENSG00000165995": "CACNB2", "ENSG00000146250": "PRSS35",
           "ENSG00000152583": "SPARCL1", "ENSG00000101347": "SAMHD1",
           "ENSG00000120129": "DUSP1", "ENSG00000189221": "MAOA",
           "ENSG00000157150": "TIMP4", "ENSG00000211445": "GPX3",
           "ENSG00000162692": "VCAM1", "ENSG00000196136": "SERPINA3"}
seen = set()
for g in res.head(10).index:
    base = g.split(".")[0]
    lab = SYMBOLS.get(base, base)
    if lab in seen:
        continue
    seen.add(lab)
    ax.annotate(lab, (res.loc[g, "log2FoldChange"], -np.log10(res.loc[g, "padj"])),
                fontsize=7, fontstyle="italic")
ax.set_xlabel("log2 fold change (Dex / untreated)")
ax.set_ylabel("-log10 adjusted p-value")
ax.set_title("airway: dexamethasone vs untreated (4 vs 4)")
plt.tight_layout(); plt.savefig(f"{OUT}/volcano_plot.png", dpi=150); plt.close()

with open(f"{OUT}/summary.txt", "w") as f:
    f.write(f"""airway (GSE52778 via recount3 SRP033351) - Dex vs untreated, design ~ cell_line + treatment
Samples: 8 (4 cell lines x 2 conditions), 4 vs 4 replicates
Genes after filter (row sum >= 10): {counts.shape[1]}
Significant FDR<{PADJ_CUTOFF}: {len(sig_fdr_only)}
Significant + |log2FC|>={LFC_CUTOFF}: {len(sig)} (up {(sig.log2FoldChange>0).sum()}, down {(sig.log2FoldChange<0).sum()})
PCA: PC1 {pca.explained_variance_ratio_[0]*100:.1f}%, PC2 {pca.explained_variance_ratio_[1]*100:.1f}%
Top gene (min padj): {SYMBOLS.get(res.index[0].split(".")[0], res.index[0])} (padj={res.padj.iloc[0]:.3g})
""")
print("\nScript 08 done.")
