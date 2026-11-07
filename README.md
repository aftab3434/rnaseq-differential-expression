# RNA-seq Differential Expression Analysis

A complete, reproducible differential expression project on real public data:
download -> QC -> DESeq2 statistics -> PCA / volcano / heatmap -> pathway enrichment.

**Dataset:** [GSE60450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450) - Fu et al. 2015,
*Nature Cell Biology* (PMID 25730472). Mouse mammary gland, luminal and basal cells,
three stages (virgin, 18.5-day pregnancy, 2-day lactation), 2 replicates each.
Raw gene-level counts are downloaded directly from GEO.

**Comparison:** lactating vs pregnant luminal cells (2 vs 2 replicates).

## Results at a glance

| Metric | Value |
|---|---|
| Genes in raw matrix | 27,179 |
| Genes after low-count filtering (>=10 counts in >=2 samples) | 15,264 |
| Genes tested (DESeq2, independent filtering) | 14,968 |
| Significant DE genes (padj < 0.05, \|log2FC\| >= 1) | **1,356** (474 up, 882 down in lactation) |
| PCA separation | PC1 61.7% + PC2 25.0% variance; replicates cluster tightly |

Top up-regulated in lactation: **Csn1s2b** (casein, a milk protein, log2FC +6.1), Mrgprf, Igfbp6.
Top down-regulated: **Cldn4** (log2FC -5.3), Lif, Gjb3.

Pathway enrichment (Enrichr, GO Biological Process + KEGG mouse) matches the known
biology of the lactation switch:
- **Up in lactation:** cholesterol / sterol / steroid biosynthesis (milk fat production)
- **Down in lactation:** mitotic sister chromatid segregation, spindle checkpoint
  (cells stop dividing and start producing milk)

| PCA | Volcano |
|---|---|
| ![PCA](results/figures/pca_plot.png) | ![Volcano](results/figures/volcano_plot.png) |

| Heatmap, top 50 DE genes | Pathway enrichment |
|---|---|
| ![Heatmap](results/figures/heatmap_top50.png) | ![Enrichment](results/figures/enrichment_dotplot.png) |

## Project structure

```
data/                  raw downloads from GEO + NCBI (script 01 fetches them)
scripts/
  01_download_data.py            dataset + gene annotation download
  02_qc.py                       library sizes, count distributions, low-count filter
  03_differential_expression.py  PyDESeq2: lactating vs pregnant luminal cells
  04_visualizations.py           PCA, volcano plot, top-50 heatmap
  05_pathway_enrichment.py       gseapy/Enrichr ORA: GO BP + KEGG mouse
results/
  figures/             the five plots the pipeline produced
  tables/              full DE results, significant genes, enrichment tables
```

## How to run

```bash
pip install -r requirements.txt
python scripts/01_download_data.py
python scripts/02_qc.py
python scripts/03_differential_expression.py
python scripts/04_visualizations.py
python scripts/05_pathway_enrichment.py   # needs internet (Enrichr API)
```

## Method notes and honest caveats

- **Statistics:** [PyDESeq2](https://pydeseq2.readthedocs.io/) (Muzellec et al. 2023,
  *Bioinformatics*), a Python re-implementation of DESeq2 matching its v1.34 defaults.
  Small numeric differences from R DESeq2 are expected and documented by its authors.
- **Replicates:** this dataset has 2 replicates per group. Benchmark studies
  (e.g. Schurch et al. 2016, *RNA*) show small replicate numbers limit statistical
  power, so this project demonstrates the analysis skill set rather than making
  strong new biological claims.
- **Enrichment:** over-representation results are hypotheses for follow-up, not proof;
  standard practice is to confirm with a second tool or literature review.
- Raw counts (never FPKM/TPM) are used for testing, and genes with very low counts
  are pre-filtered, per DESeq2/edgeR standard practice.

## License

MIT (see LICENSE). Data: public GEO dataset, please cite Fu et al. 2015 if reused.
