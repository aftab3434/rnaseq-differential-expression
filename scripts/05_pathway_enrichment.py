"""Step 05 - Pathway enrichment of the significant genes (gseapy / Enrichr).

Over-representation analysis of genes up- and down-regulated in lactation,
against GO Biological Process and KEGG (mouse) via the Enrichr API.
Enrichment results are hypotheses for follow-up, not proof.

Outputs:
  results/tables/enrichment_up_in_lactation.csv
  results/tables/enrichment_down_in_lactation.csv
  results/figures/enrichment_dotplot.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import gseapy as gp

ROOT = Path(__file__).resolve().parent.parent
FIGS, TABLES = ROOT / "results" / "figures", ROOT / "results" / "tables"

sig = pd.read_csv(TABLES / "de_genes_significant.csv", index_col=0)
sig = sig[sig["symbol"].notna() & (sig["symbol"] != "")]
up = sig.loc[sig["log2FoldChange"] > 0, "symbol"].tolist()
down = sig.loc[sig["log2FoldChange"] < 0, "symbol"].tolist()
print(f"symbols: {len(up)} up, {len(down)} down")

def enrich(genes, label):
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=["GO_Biological_Process_2025", "KEGG_2019_Mouse"],
        organism="mouse",
        outdir=None,
    )
    df = enr.results.sort_values("Adjusted P-value")
    df.to_csv(TABLES / f"enrichment_{label}.csv", index=False)
    return df

up_res = enrich(up, "up_in_lactation")
down_res = enrich(down, "down_in_lactation")

# dotplot: top 8 terms per direction by adjusted p-value
def top(df, n=8):
    t = df.nsmallest(n, "Adjusted P-value").copy()
    t["Term"] = t["Term"].str.replace(r"\s*\(.*\)$", "", regex=True)
    return t

fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=False)
for ax, df, title, color in [
    (axes[0], top(up_res), "Up in lactation", "#d62728"),
    (axes[1], top(down_res), "Down in lactation", "#1f77b4"),
]:
    df = df.iloc[::-1]
    ax.scatter(-df["Adjusted P-value"].map(lambda p: __import__("numpy").log10(p)),
               range(len(df)), s=df["Odds Ratio"].clip(upper=15) * 12, color=color, alpha=0.75)
    ax.set_yticks(range(len(df)))
    import textwrap
    ax.set_yticklabels([textwrap.fill(t, 38) for t in df["Term"]], fontsize=8)
    ax.set_xlabel("-log10 adjusted p-value")
    ax.set_title(title)
fig.suptitle("Pathway enrichment (Enrichr: GO Biological Process + KEGG mouse)")
fig.tight_layout()
fig.savefig(FIGS / "enrichment_dotplot.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("\ntop enriched terms, up in lactation:")
print(up_res.head(5)[["Gene_set", "Term", "Adjusted P-value"]].to_string(index=False))
print("\ntop enriched terms, down in lactation:")
print(down_res.head(5)[["Gene_set", "Term", "Adjusted P-value"]].to_string(index=False))
