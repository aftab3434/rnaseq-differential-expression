"""Step 04 - Visualization: PCA, volcano plot, heatmap of top genes.

Reads the DESeq2-normalized counts and the DE results from steps 02-03.
Writes three figures into results/figures/.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parent.parent
FIGS = ROOT / "results" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

norm = pd.read_csv(ROOT / "data" / "processed" / "normed_counts_luminal.csv", index_col=0)
meta = pd.read_csv(ROOT / "data" / "processed" / "sample_metadata.csv")
meta = meta[meta["cell_type"] == "luminal"].set_index("sample").loc[norm.columns]
res = pd.read_csv(ROOT / "results" / "tables" / "de_results_lactating_vs_pregnant_full.csv", index_col=0)

# ---------- PCA (log2 normalized counts, numpy SVD - no extra dependencies) ----------
logc = np.log2(norm + 1)
X = (logc - logc.mean(axis=1).values[:, None]).T.values  # samples x genes, gene-centred
U, S, _ = np.linalg.svd(X, full_matrices=False)
var = S**2 / (S**2).sum() * 100

fig, ax = plt.subplots(figsize=(6, 5))
for stage, color in [("virgin", "#1f77b4"), ("pregnant", "#2ca02c"), ("lactating", "#d62728")]:
    m = (meta["stage"] == stage).values
    ax.scatter(U[m, 0] * S[0], U[m, 1] * S[1], s=90, label=stage, color=color)
    for x, y, name in zip(U[m, 0] * S[0], U[m, 1] * S[1], meta.index[m]):
        ax.annotate(name.replace("MCL1-", ""), (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8)
ax.set_xlabel(f"PC1 ({var[0]:.1f}% variance)")
ax.set_ylabel(f"PC2 ({var[1]:.1f}% variance)")
ax.set_title("PCA of luminal samples (log2 normalized counts)")
ax.legend()
fig.tight_layout()
fig.savefig(FIGS / "pca_plot.png", dpi=150)
plt.close(fig)
print(f"PCA: PC1 {var[0]:.1f}%, PC2 {var[1]:.1f}%")

# ---------- Volcano plot ----------
res = res.dropna(subset=["padj"])
res["neglog10padj"] = -np.log10(res["padj"])
sig = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)

fig, ax = plt.subplots(figsize=(7, 5.5))
ax.scatter(res.loc[~sig, "log2FoldChange"], res.loc[~sig, "neglog10padj"], s=4, color="grey", alpha=0.4, label="not significant")
ax.scatter(res.loc[sig & (res["log2FoldChange"] > 0), "log2FoldChange"], res.loc[sig & (res["log2FoldChange"] > 0), "neglog10padj"],
           s=6, color="#d62728", alpha=0.7, label="up in lactation")
ax.scatter(res.loc[sig & (res["log2FoldChange"] < 0), "log2FoldChange"], res.loc[sig & (res["log2FoldChange"] < 0), "neglog10padj"],
           s=6, color="#1f77b4", alpha=0.7, label="down in lactation")
for _, row in res.nsmallest(10, "padj").iterrows():
    if row["symbol"]:
        ax.annotate(row["symbol"], (row["log2FoldChange"], row["neglog10padj"]),
                    textcoords="offset points", xytext=(4, 4), fontsize=8)
ax.axvline(1, ls="--", lw=0.6, color="black"); ax.axvline(-1, ls="--", lw=0.6, color="black")
ax.axhline(-np.log10(0.05), ls="--", lw=0.6, color="black")
ax.set_xlabel("log2 fold change (lactating vs pregnant)")
ax.set_ylabel("-log10 adjusted p-value")
ax.set_title("Volcano plot: lactating vs pregnant luminal cells")
ax.legend(markerscale=2)
fig.tight_layout()
fig.savefig(FIGS / "volcano_plot.png", dpi=150)
plt.close(fig)

# ---------- Heatmap of top 50 DE genes ----------
top50 = res.nsmallest(50, "padj")
mat = np.log2(norm.loc[top50.index] + 1)
z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
z.index = [s if s else str(g) for s, g in zip(top50["symbol"], top50.index)]

col_colors = meta["stage"].map({"virgin": "#1f77b4", "pregnant": "#2ca02c", "lactating": "#d62728"})
g = sns.clustermap(z, cmap="RdBu_r", center=0, col_cluster=False,
                   col_colors=col_colors, figsize=(8, 10),
                   yticklabels=True, cbar_kws={"label": "row z-score"})
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=8)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)
g.fig.suptitle("Top 50 DE genes (by adjusted p-value), lactating vs pregnant", y=1.01)
g.savefig(FIGS / "heatmap_top50.png", dpi=150, bbox_inches="tight")
plt.close("all")
print("wrote pca_plot.png, volcano_plot.png, heatmap_top50.png")
