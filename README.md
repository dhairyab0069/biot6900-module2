# Multi-Omics Target Prioritization in Alzheimer's Disease

**BIOT 6900 · Module 2 (Assignment 2): Independent Multi-Omics Target Discovery**
Dhairya Bhatia · Northeastern University

This project ranks candidate drug targets for **Alzheimer's disease (AD)** by combining three
layers of public data: bulk brain transcriptomics, brain proteomics and GWAS. Each layer is
reduced to one row per gene and joined on the gene symbol. RNA and protein are checked for
direction agreement, and genes are ranked by an equal-weighted multi-evidence score.

## Repository contents

| Path | Description |
|---|---|
| `BIOT6900_Module2_Starter.ipynb` | Analysis notebook (Parts 1 to 3). Runs top to bottom with `Restart & Run All`. |
| `scripts/get_data.py` | Downloads every dataset and writes the processed tables to `data/`. |
| `data/` | Processed input tables (see [Data sources](#data-sources)). |
| `data/targets_ad.csv` | **Deliverable.** Top 15 ranked AD targets. |
| `data/targets_ad_full.csv` | Full ranked list (1,115 genes), passed on to Week 3. |
| `REPORT.md` / `REPORT.pdf` | Written report. |

## Reproducing the analysis

Tested with Python 3.13, pandas 3.0, NumPy 2.5 and SciPy 1.18.

```bash
python3 scripts/get_data.py        # ~140 MB download on first run, cached in .cache/
jupyter lab BIOT6900_Module2_Starter.ipynb
# then: Kernel > Restart & Run All
```

## Data sources

All data was sourced independently (the Canvas matrices were not used). Every dataset is
publicly available.

| Layer | Dataset | Source | Access | Genes |
|---|---|---|---|---|
| Transcriptomics | AMP-AD harmonized RNA-seq, AD vs. control. ROSMAP, Mayo and MSBB cohorts, 9 brain regions. | [Agora](https://agora.adknowledgeportal.org) ([API](https://agora.adknowledgeportal.org/api/v1/genes/comparison?category=RNA%20-%20Differential%20Expression&subCategory=AD%20Diagnosis%20(males%20and%20females))) | Open | 20,539 |
| Proteomics | AMP-AD TMT mass-spec proteomics, AD vs. control. Dorsolateral prefrontal cortex. | [Agora](https://agora.adknowledgeportal.org) ([API](https://agora.adknowledgeportal.org/api/v1/genes/comparison?category=Protein%20-%20Differential%20Expression&subCategory=TMT)) | Open | 8,254 |
| Genomics | All GWAS associations mapped to Alzheimer disease (`MONDO_0004975`). | [GWAS Catalog](https://www.ebi.ac.uk/gwas/efotraits/MONDO_0004975) ([release](https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/)) | Open | 3,558 |
| Part 1 demo | CPTAC breast cancer RNA, protein and somatic mutation. | [LinkedOmics CPTAC-BRCA](https://www.linkedomics.org/data_download/CPTAC-BRCA/) | Open | 23,121 / 12,621 / 9,448 |

**Access.** Agora publishes gene-level summary statistics openly. The individual-level AMP-AD
data behind them is hosted on Synapse and requires a Data Use Certificate (application). It was
not needed and was not used.

**Study design.** The layers are **unmatched**. Each is a per-gene summary with no sample
identifiers, so integration is at the gene level rather than the sample level.

## Methods

1. **Preprocessing** (`scripts/get_data.py`)
   - RNA: for each gene, keep the brain region with the largest absolute log2 fold change.
     Cerebellum is excluded because it is relatively spared in AD.
   - Protein: collapse multiple UniProt isoforms per gene symbol to the largest effect.
   - GWAS: split multi-gene associations into one row per gene and keep each gene's strongest
     association (maximum −log10 p).
   - RNA and protein p-values are Agora's multiple-testing-adjusted values (`adj_p_val`).
2. **Harmonization.** Inner join on the HGNC gene symbol.
3. **Concordance.** `concordant = sign(RNA lfc) == sign(protein lfc)`.
4. **Scoring.** Per-layer magnitude (|RNA lfc|, |protein lfc|, −log10 p), rank-percentile
   normalized, then combined with equal weights (1/3 each).
5. **Export.** Sort by score and write the top 15 and the full list to `data/`.

## Key results

- **1,115 genes** are present in all three layers. The proteomics layer limits coverage.
- **635 / 1,115 (57%)** are sign-concordant between RNA and protein.
- Known AD risk genes recovered: **APOE (rank 11)**, HLA-DRB1 (2), CLU (31), BIN1 (53),
  PICALM (83), SORL1 (220), ABCA7 (585).
- Top hits cluster into microglial activation (HLA-DRA, HLA-DRB1, ITGAX), lipid transport
  (APOC1, APOE, APOB), synaptic loss (RPH3A, NRXN3, SHANK2) and vascular genes (PECAM1, APLNR).

See `REPORT.pdf` for interpretation.

## Limitations

- **No per-patient claims.** Unmatched summary data cannot show RNA/protein coupling within
  individuals, unlike the sample-matched CPTAC analysis in Part 1.
- **Proteomic coverage.** TREM2, CD33 and MS4A6A have strong RNA and GWAS evidence but are not
  quantified in the TMT data, so the inner join drops them.
- **Region mismatch.** RNA uses the strongest region per gene; protein is DLPFC only. Restricting
  RNA to DLPFC weakens known signal (APOE drops from rank 11 to 110).
- **Linkage at 19q13.** APOC1 and TOMM40 inherit much of their GWAS signal from APOE.

## Issues encountered

- The GWAS Catalog download API (`api/search/downloads`) returned empty files or HTTP 500 errors.
  The full association release is downloaded from the FTP site and filtered locally instead.
- The trait ID `EFO_0000249` returns no rows; the Catalog now maps Alzheimer disease to
  `MONDO_0004975`.
- LinkedOmics returns HTTP 403 to Python's default client, so the script sends a browser
  `User-Agent` header.
- The CPTAC-BRCA download has no tumor/normal labels, so Part 1 ranks by overall abundance
  (TP53 is rank 108 of 6,306).
