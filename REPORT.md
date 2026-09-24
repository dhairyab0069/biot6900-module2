# BIOT 6900 · Module 2: Finding Drug Targets for Alzheimer's Disease

**Dhairya Bhatia** · bhatia.dh@northeastern.edu

## Key terms

- **RNA:** how much a gene is being "read" in a tissue.
- **Protein:** how much of the gene's final product is present.
- **GWAS:** studies that find DNA variants linked to a disease.
- **Fold change (lfc):** how much a gene changes in AD vs. healthy brains. Above 0 = up, below 0 = down.
- **−log10 p:** how strong a genetic link is. Bigger = stronger.
- **Concordant:** RNA and protein move in the same direction.

## 1. Data

I found and downloaded all the data myself with one script (`scripts/get_data.py`). All sources
are free, with no login or application needed.

| Layer | Source | Access | Genes |
|-------------------|-------------------------------------|--------|--------|
| RNA | Agora (agora.adknowledgeportal.org) | Open | 20,539 |
| Protein | Agora (agora.adknowledgeportal.org) | Open | 8,254 |
| Genetics | GWAS Catalog (ebi.ac.uk/gwas) | Open | 3,558 |

- **RNA:** AD vs. healthy brains, from 3 donor groups (ROSMAP, Mayo, MSBB) and 9 brain regions.
- **Protein:** AD vs. healthy brains, from the prefrontal cortex.
- **Genetics:** every published DNA variant linked to Alzheimer's.

The data is unmatched. Each file gives one number per gene, not per person, so I joined the
layers by gene name.

## 2. Steps

1. Load the three files.
2. For RNA, keep each gene's biggest change across brain regions (skipping the cerebellum, which
   Alzheimer's mostly spares).
3. Keep only genes found in all three layers. **1,115 genes** made it.
4. Mark a gene **concordant** if RNA and protein move the same way. **635 of 1,115 (57%)** do.
5. Rank each layer from 0 to 1 and add the three ranks with equal weight (1/3 each).
6. Sort by score. Save the top 15 to `data/targets_ad.csv` and all 1,115 to
   `data/targets_ad_full.csv`.

## 3. Results

### Top 15 genes

| Rank | Gene | RNA | Protein | Genetics | Concordant | Score |
|---|---|---|---|---|---|---|
| 1 | HLA-DRA | 0.80 | 0.14 | 16.7 | Yes | 0.957 |
| 2 | HLA-DRB1 | 0.66 | 0.10 | 22.5 | Yes | 0.955 |
| 3 | APOC1 | 0.51 | −0.23 | 672.7 | No | 0.950 |
| 4 | PECAM1 | 0.86 | 0.09 | 16.5 | Yes | 0.939 |
| 5 | MTSS2 | 0.64 | −0.10 | 13.5 | No | 0.913 |
| 6 | ITGAX | 0.70 | 0.12 | 11.7 | Yes | 0.912 |
| 7 | GMPR | 1.04 | 0.23 | 10.1 | Yes | 0.910 |
| 8 | RPH3A | −1.16 | −0.13 | 10.0 | Yes | 0.905 |
| 9 | NRXN3 | −0.64 | 0.06 | 22.0 | No | 0.904 |
| 10 | APLNR | 1.28 | 0.10 | 10.5 | Yes | 0.902 |
| 11 | **APOE** | 0.38 | 0.11 | 496.4 | Yes | 0.902 |
| 12 | SHANK2 | −0.55 | −0.06 | 30.0 | Yes | 0.899 |
| 13 | APOB | 0.40 | 0.13 | 19.0 | Yes | 0.897 |
| 14 | MOG | 0.47 | −0.09 | 17.3 | No | 0.896 |
| 15 | PLEC | 0.60 | 0.26 | 10.5 | Yes | 0.895 |

*RNA and Protein = fold change. Genetics = −log10 p.*

## 4. What the results mean

### Why Alzheimer's and why this data

- Alzheimer's is the most common cause of dementia and has very few treatments.
- It follows on from Module 1, which used the APOE variant `rs7412`.
- Agora combines 3 donor groups and 9 brain regions, so it is bigger than a single GEO study.
- The GWAS Catalog collects every published Alzheimer's genetics result.

### Do the top genes make sense?

Yes. Known Alzheimer's risk genes show up:

| Gene | Rank | Note |
|----------|------|--------------------------------------------------------------------|
| APOE | 11 | Biggest known genetic risk for Alzheimer's. Up at RNA and protein. |
| HLA-DRB1 | 2 | Known risk gene. Immune system. |
| CLU | 31 | Known risk gene. Fat transport. |
| BIN1 | 53 | Known risk gene. |
| PICALM | 83 | Known risk gene. |
| SORL1 | 220 | Known risk gene. |

The top 15 genes fall into four groups that are all known to change in Alzheimer's:

1. **Brain immune cells:** HLA-DRA, HLA-DRB1, ITGAX (all up).
2. **Fat transport:** APOC1, APOE, APOB.
3. **Synapses (nerve connections):** RPH3A, NRXN3, SHANK2 (all down, which fits synapse loss).
4. **Blood vessels:** PECAM1, APLNR.

#### Surprise hit: GMPR (rank 7)

GMPR is not a known Alzheimer's gene, but it goes up at both RNA
(+1.04) and protein (+0.23), and both results are strong. That makes it worth a closer look.

#### Missing genes

TREM2, CD33 and MS4A6A are known Alzheimer's genes but are not in the list.
The protein study did not measure them, so they were dropped in step 3.

#### A warning about APOC1

APOC1 (rank 3) sits right next to APOE on chromosome 19. Its huge genetic
score (673) is probably mostly APOE's signal, not its own.

### A gene where RNA and protein disagree: APOC1

| | Change | Brain region |
|---|---|---|
| RNA | up (+0.51) | parahippocampal gyrus |
| Protein | down (−0.23) | prefrontal cortex |

Both changes are statistically strong. Possible reasons:

1. RNA and protein were measured in **different brain regions**.
2. APOC1 is a protein that cells **send out**, so less stays in the tissue.
3. The protein may be **broken down faster** or **trapped in plaques**.
4. There may be **more immune cells** in the tissue, raising RNA without each cell making more.

Not every "No" is real. BIN1 is marked not concordant, but its protein change is tiny and not
significant (p = 0.43), so that is probably just noise.

### What we can and cannot say

#### What we can say

- These genes have several separate kinds of evidence (RNA, protein and genetics).
- The three layers barely overlap, so a high score is not one signal counted three times.

#### What we cannot say

1. We can't say anything per patient. In Part 1 (breast cancer), RNA and protein came from the same
   tumors, so we could compare them person by person. Here the data is unmatched.
2. We can't say what causes what. A genetic link points to a region of DNA, not a proven cause.
3. We can't say the missing genes don't matter. TREM2 was dropped only because protein wasn't measured.
4. We can't give exact effect sizes. Taking the biggest RNA change per gene makes RNA effects look larger.
   Only 15% of genes have a significant protein change.

### Weighting

I kept equal weights (1/3 each). Ranking each layer first stops APOE's huge genetic score from
taking over. Protein could get less weight because its changes are small, but I kept equal
weights to match the course default.

### Conclusion

Combining RNA, protein and genetics found known Alzheimer's genes like APOE, HLA-DRB1 and CLU
near the top. That suggests the method works. The ranked list is ready for Week 3, and genes with
strong evidence in all three layers, such as APOE, the HLA genes and GMPR, are the best places to
start.

If I did this again, I would:

1. Not require all three layers, so genes like TREM2 are not lost just because protein wasn't
   measured.
2. Use RNA and protein from the same brain region.
3. Count a gene as concordant only when both changes are significant.
4. Treat APOE and its neighbours (APOC1, TOMM40) as one genetic signal.
