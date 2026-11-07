"""Step 01 - Download the public dataset (GSE60450).

Downloads three files into data/:
  1. GSE60450_Lactation-GenewiseCounts.txt.gz - gene-level raw counts (GEO supplementary)
  2. GSE60450_series_matrix.txt.gz            - sample metadata (GEO series matrix)
  3. Mus_musculus.gene_info.gz                - NCBI GeneID -> symbol mapping

Study: Fu et al. 2015, Nature Cell Biology (PMID 25730472).
Mouse mammary gland, luminal and basal cells, three stages
(virgin, 18.5-day pregnancy, 2-day lactation), 2 replicates each.
"""
import gzip
import shutil
import urllib.request
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
DATA.mkdir(exist_ok=True)

FILES = {
    "GSE60450_Lactation-GenewiseCounts.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz",
    "GSE60450_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/matrix/GSE60450_series_matrix.txt.gz",
    "Mus_musculus.gene_info.gz":
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
}

for name, url in FILES.items():
    dest = DATA / name
    if dest.exists():
        print(f"already present: {name}")
        continue
    print(f"downloading {name} ...")
    urllib.request.urlretrieve(url, dest)
    print(f"  saved {dest} ({dest.stat().st_size/1e6:.1f} MB)")

# sanity check: the count file must open and have the expected header
with gzip.open(DATA / "GSE60450_Lactation-GenewiseCounts.txt.gz", "rt") as fh:
    header = fh.readline().strip().split("\t")
assert header[0] == "EntrezGeneID", f"unexpected count file header: {header[:3]}"
print(f"count file OK: {len(header) - 2} samples + EntrezGeneID + Length")
