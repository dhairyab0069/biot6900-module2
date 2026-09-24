import io
import json
import os
import urllib.request
import zipfile

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")
CACHE = os.path.join(REPO, ".cache")

LINKEDOMICS = "https://www.linkedomics.org/data_download/CPTAC-BRCA/"
AGORA = "https://agora.adknowledgeportal.org/api/v1/genes/comparison"
AGORA_RNA = AGORA + "?category=RNA%20-%20Differential%20Expression&subCategory=AD%20Diagnosis%20(males%20and%20females)"
AGORA_PROTEIN = AGORA + "?category=Protein%20-%20Differential%20Expression&subCategory=TMT"
GWAS_CATALOG = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
                "gwas-catalog-associations_ontology-annotated-full.zip")

# The GWAS Catalog's ID for Alzheimer's disease. (The older ID, EFO_0000249, now finds nothing.)
ALZHEIMERS_ID = "MONDO_0004975"

# Brain regions to skip. The cerebellum (CBE) is mostly spared by Alzheimer's.
SKIP_REGIONS = {"CBE"}


def download(url, filename):
    """Download a file into .cache/, or reuse it if it is already there."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, filename)
    if os.path.exists(path):
        print(f"  already downloaded: {filename}")
        return path
    print(f"  downloading {filename} ...")
    # LinkedOmics blocks requests that don't look like they come from a web browser.
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=600) as response:
        content = response.read()
    with open(path + ".part", "wb") as f:   # rename only once fully downloaded
        f.write(content)
    os.rename(path + ".part", path)
    print(f"  saved {os.path.getsize(path) / 1e6:.1f} MB")
    return path


def save(table, filename, keep_row_names=False):
    table.to_csv(os.path.join(DATA, filename), sep="\t", index=keep_row_names)
    print(f"  wrote data/{filename}: {len(table)} genes")


# ---------------------------------------------------------------- breast cancer

def read_linkedomics(filename):
    """Read a LinkedOmics table: one row per gene, one column per tumor."""
    return pd.read_csv(download(LINKEDOMICS + filename, filename), sep="\t", index_col=0)


def get_breast_cancer():
    print("\nBreast cancer (LinkedOmics)")
    save(read_linkedomics("HS_CPTAC_BRCA_2018_RNA_GENE.cct"),
         "cptac_brca_rna.tsv", keep_row_names=True)
    save(read_linkedomics("HS_CPTAC_BRCA_2018_Proteome_Ratio_Norm_gene_Median.cct"),
         "cptac_brca_protein.tsv", keep_row_names=True)

    # The mutation table is 1 (mutated) or 0 (not) per tumor. Averaging a row gives the
    # fraction of tumors where that gene is mutated (rounded to 6 decimal places).
    mutated = read_linkedomics("HS_CPTAC_BRCA_2018_MUT_GENE.cbt")
    save(mutated.mean(axis=1).round(6).rename("mut_freq").to_frame(),
         "cptac_brca_mutation.tsv", keep_row_names=True)


# ------------------------------------------------------------------ Alzheimer's

def read_agora(url, filename):
    """Agora lists each gene with a result for every brain region it was measured in.
    Turn that into a simple table with one row per gene per region."""
    with open(download(url, filename)) as f:
        genes = json.load(f)["items"]
    rows = [{"gene": g["hgnc_symbol"], "tissue": region["name"],
             "log2fc": region["logfc"], "pval": region["adj_p_val"]}
            for g in genes for region in g["tissues"]]
    table = pd.DataFrame(rows).dropna(subset=["gene", "log2fc"])
    return table[table["gene"].str.strip() != ""]   # a few entries have no gene name


def biggest_change_per_gene(table):
    """Some genes appear several times (several brain regions, or several versions of the
    same protein). Keep only the row with the biggest change, up or down."""
    best_rows = table.groupby("gene")["log2fc"].apply(lambda fc: fc.abs().idxmax())
    return table.loc[best_rows].sort_values("gene")


def get_alzheimers_rna():
    print("\nAlzheimer's RNA (Agora)")
    table = read_agora(AGORA_RNA, "agora_rna_de.json")
    table = table[~table["tissue"].isin(SKIP_REGIONS)]
    save(biggest_change_per_gene(table)[["gene", "log2fc", "pval", "tissue"]],
         "ad_transcriptomics.tsv")


def get_alzheimers_protein():
    print("\nAlzheimer's protein (Agora)")
    table = read_agora(AGORA_PROTEIN, "agora_tmt_de.json")
    save(biggest_change_per_gene(table)[["gene", "log2fc", "pval"]], "ad_proteomics.tsv")


def get_alzheimers_genetics():
    print("\nAlzheimer's genetics (GWAS Catalog)")
    path = download(GWAS_CATALOG, "gwas_catalog_associations.zip")

    # The Catalog covers every disease and is large, so read it in pieces and keep only
    # the Alzheimer's rows.
    columns = ["MAPPED_TRAIT_URI", "MAPPED_GENE", "PVALUE_MLOG"]
    pieces = []
    with zipfile.ZipFile(path) as z:
        inner = next(name for name in z.namelist() if name.endswith(".tsv"))
        with z.open(inner) as f:
            text = io.TextIOWrapper(f, encoding="utf-8", errors="replace")
            for piece in pd.read_csv(text, sep="\t", usecols=columns,
                                     low_memory=False, chunksize=200_000):
                is_ad = piece["MAPPED_TRAIT_URI"].astype(str).str.contains(ALZHEIMERS_ID)
                pieces.append(piece[is_ad])
    ad = pd.concat(pieces).dropna(subset=["MAPPED_GENE", "PVALUE_MLOG"])
    print(f"  {len(ad)} Alzheimer's results found")

    # One result can name several genes, e.g. "APOE, APOC1" or "BIN1 - CYP27C1".
    # Split those into one row per gene. (Don't split on a plain "-": it's part of real
    # gene names like HLA-DRA.)
    ad["gene"] = ad["MAPPED_GENE"].astype(str).str.split(r"\s*[,;]\s*|\s+-\s+", regex=True)
    ad = ad.explode("gene")
    ad["gene"] = ad["gene"].str.strip()
    ad = ad[ad["gene"].str.match(r"^[A-Za-z][A-Za-z0-9._-]*$", na=False)]

    # PVALUE_MLOG is -log10(p): bigger = stronger link. Keep each gene's strongest result.
    strongest = ad.groupby("gene", as_index=False)["PVALUE_MLOG"].max()
    save(strongest.rename(columns={"PVALUE_MLOG": "neglog10p"}), "ad_gwas.tsv")


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    get_breast_cancer()
    get_alzheimers_rna()
    get_alzheimers_protein()
    get_alzheimers_genetics()
    print("\nDone. All six files are in data/.")
