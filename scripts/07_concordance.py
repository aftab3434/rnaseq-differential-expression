"""Step 07 - Concordance between the two independent cell types.

If the lactation signature is real, genes should change in the same
direction in luminal and basal cells. Plots luminal vs basal log2FC and
reports the Spearman correlation, genome-wide and for significant genes.

Outputs: results/figures/concordance_scatter.png and printed statistics.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
TABLES, FIGS = ROOT / "results" / "tables", ROOT / "results" / "figures"

lum = pd.read_csv(TABLES / "de_results_lactating_vs_pregnant_full.csv", index_col=0)
bas = pd.read_csv(TABLES / "de_results_basal_lactating_vs_pregnant_full.csv", index_col=0)

m = lum[["symbol", "log2FoldChange", "padj"]].join(
    bas[["log2FoldChange", "padj"]], lsuffix="_luminal", rsuffix="_basal", how="inner"
).dropna()
print(f"genes tested in both contrasts: {len(m):,}")

rho_all, p_all = stats.spearmanr(m["log2FoldChange_luminal"], m["log2FoldChange_basal"])
sig_either = m[(m["padj_luminal"] < 0.05) | (m["padj_basal"] < 0.05)]
rho_sig, _ = stats.spearmanr(sig_either["log2FoldChange_luminal"], sig_either["log2FoldChange_basal"])
same_dir = ((sig_either["log2FoldChange_luminal"] * sig_either["log2FoldChange_basal"]) > 0).mean()
print(f"Spearman correlation of log2FC, all genes: rho = {rho_all:.3f}")
print(f"Spearman correlation, genes significant in either contrast: rho = {rho_sig:.3f} (n={len(sig_either):,})")
print(f"same-direction fraction among those: {same_dir:.1%}")

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(m["log2FoldChange_luminal"], m["log2FoldChange_basal"], s=3, color="grey", alpha=0.3, label="all tested genes")
ax.scatter(sig_either["log2FoldChange_luminal"], sig_either["log2FoldChange_basal"],
           s=6, color="#d62728", alpha=0.6, label="significant in either contrast")
ax.axhline(0, lw=0.6, color="black"); ax.axvline(0, lw=0.6, color="black")
ax.set_xlabel("log2FC in luminal cells (lactating vs pregnant)")
ax.set_ylabel("log2FC in basal cells (lactating vs pregnant)")
ax.set_title(f"Lactation signature replicates across cell types\nSpearman rho = {rho_all:.2f} (all genes)")
ax.legend(markerscale=2)
fig.tight_layout()
fig.savefig(FIGS / "concordance_scatter.png", dpi=150)
plt.close(fig)
print("wrote concordance_scatter.png")
