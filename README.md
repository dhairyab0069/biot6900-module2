# BIOT 6900 · Module 2 — Multi-Omics Target Identification & Validation

Dhairya Bhatia

**Disease: Alzheimer's disease (AD)**

## Deliverables

- `BIOT6900_Module2_Starter.ipynb` — runs top to bottom (`Restart & Run All`, verified).
- `data/targets_ad.csv` — the top 15 ranked Alzheimer's targets.
- `data/targets_ad_full.csv` — the full ranked list (1,115 genes), for Week 3.
- `REPORT.md` — the written report.
- `scripts/get_data.py` — downloads every dataset below and saves it into `data/`.

## Datasets

All data was found and downloaded independently (not the Canvas files). Every source below is
free and needs no login or application.

| Layer | File | Dataset | Source link | Access |
|---|---|---|---|---|
| Transcriptomics | `data/ad_transcriptomics.tsv` | AMP-AD harmonized RNA-seq, AD vs. control (ROSMAP, Mayo, MSBB cohorts; 9 brain regions) | [Agora](https://agora.adknowledgeportal.org) · [API](https://agora.adknowledgeportal.org/api/v1/genes/comparison?category=RNA%20-%20Differential%20Expression&subCategory=AD%20Diagnosis%20(males%20and%20females)) | Open |
| Proteomics | `data/ad_proteomics.tsv` | AMP-AD TMT mass-spec proteomics, AD vs. control (DLPFC) | [Agora](https://agora.adknowledgeportal.org) · [API](https://agora.adknowledgeportal.org/api/v1/genes/comparison?category=Protein%20-%20Differential%20Expression&subCategory=TMT) | Open |
| Genomics | `data/ad_gwas.tsv` | All GWAS associations mapped to Alzheimer disease (`MONDO_0004975`) | [GWAS Catalog trait page](https://www.ebi.ac.uk/gwas/efotraits/MONDO_0004975) · [full release](https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/) | Open |
| Breast cancer (Part 1 demo) | `data/cptac_brca_*.tsv` | CPTAC breast cancer RNA, protein, mutation | [LinkedOmics CPTAC-BRCA](https://www.linkedomics.org/data_download/CPTAC-BRCA/) | Open |

**Access note.** Agora serves *gene-level summary statistics* (fold change and p-value per gene),
which are open. The *individual-patient* AMP-AD data behind them lives on Synapse and requires a
signed data use agreement — that is not needed here and was not used.

**Matched or unmatched?** Unmatched. Each layer is a per-gene summary (one fold change or
p-value per gene), with no sample IDs, so RNA and protein can't be paired person by person.
The GWAS results come from separate genetic studies. The RNA and protein cohorts may share some
donors (both draw on AMP-AD brain banks), but that can't be used here. So the layers are
integrated at the gene level, not the sample level.

**p-values.** Agora's `pval` for RNA and protein is the multiple-testing-adjusted p-value
(`adj_p_val`).

Two judgment calls worth flagging:

- **RNA brain region.** Agora reports DE per brain region. I take, per gene, the
  largest-magnitude effect across the eight AD-affected regions and drop cerebellum, which is
  relatively spared in AD. Restricting instead to DLPFC alone (to tissue-match the proteomics)
  gives a similar join size but much weaker signal — APOE falls from rank 11 to 110.
- **Gene symbol collapsing.** A TMT symbol can carry several UniProt isoforms, and one GWAS
  association can map to several genes. Both are collapsed to one row per gene by strongest
  effect.

To regenerate everything: `python3 scripts/get_data.py` (downloads ~140 MB once, kept in `.cache/`).

## Results

1,115 genes survive the 3-way join — the proteome is the bottleneck, since mass spec quantifies
far fewer genes than RNA-seq measures. 635 of them (57%) are sign-concordant between RNA and
protein.

Canonical AD genes recovered: **APOE 11**, CLU 31, BIN1 53, PICALM 83, SORL1 220, ABCA7 585.

**TREM2 is absent from the ranked list.** It clears the RNA and GWAS layers comfortably but is
not quantified in the TMT proteome — a low-abundance microglial membrane protein — so the inner
join drops it. That is a coverage limitation of the proteomic layer, not a claim that TREM2 is
uninteresting, and it is the clearest argument in this dataset against requiring all three
layers before a gene can be nominated.

## What didn't work

- The GWAS Catalog's `api/search/downloads` endpoint returns a header-only file or a Tomcat 500
  error, so the script pulls the full association release from the FTP site and filters locally.
- Filtering that release on `EFO_0000249` (the trait ID named in the notebook's APOE anchor
  context) returns **zero** rows; the Catalog has remapped Alzheimer disease to
  `MONDO_0004975`.
- Part 2's real CPTAC download has no tumor/normal labels, so Part 1 ranks on overall abundance
  rather than a tumor-vs-normal effect. TP53 lands at rank 108 of 6,306 rather than near the
  top — expected, and noted in the Part 2 cell.
- LinkedOmics rejects Python's default web client (HTTP 403) even though a browser or `curl`
  works, so the script sends a browser-style `User-Agent` header.
