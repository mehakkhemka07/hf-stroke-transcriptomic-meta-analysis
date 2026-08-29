# Consensus Signature Markers Interconnecting Heart Failure and Ischemic Stroke

A multi-dataset transcriptomic meta-analysis identifying shared gene expression signatures between heart failure (HF) and ischemic stroke (IS) across 7 GEO datasets, using DESeq2, a Priority Index (PI) scoring system, K-means clustering, GO/KEGG enrichment, and STRING PPI network analysis to prioritize cross-disease hub genes.

## Repository Structure

```
hf-stroke-transcriptomic-meta-analysis/
├── 01_datasets/
│   ├── download_metadata.R            # GEOquery script — pulls sample metadata for all 7 GEO datasets
│   └── data_sources.md                # GEO accessions + raw count supplementary file names used
├── data/
│   ├── GSE122709_all.counts
│   ├── GSE140275_raw_counts_GRCh38.p13_NCBI.tsv
│   ├── GSE56267_raw_counts_GRCh38.p13_NCBI.tsv
│   ├── GSE178764_raw_counts_GRCh38.p13_NCBI.tsv
│   ├── GSE116250_raw_counts_GRCh38.p13_NCBI.tsv
│   ├── GSE120852_raw_counts_GRCh38.p13_NCBI.tsv
│   ├── GSE135055_raw_counts_GRCh38.p13_NCBI.tsv
│   ├── GSE122709_sample_metadata.csv
│   ├── GSE140275_sample_metadata.csv
│   ├── GSE56267_sample_metadata.csv
│   ├── GSE178764_sample_metadata.csv
│   ├── GSE116250_sample_metadata.csv
│   ├── GSE120852_sample_metadata.csv
│   ├── GSE135055_sample_metadata.csv
│   └── gene_annotation.csv            # GeneID → Symbol lookup used in run_deseq2.R
├── 02_deseq2/
│   ├── scripts/
│   │   └── run_deseq2.R               # Runs DESeq2 across all 7 datasets, split by disease group
│   └── results/
│       ├── GSE122709_deseq2.csv
│       ├── GSE140275_deseq2.csv
│       ├── GSE56267_deseq2.csv
│       ├── GSE178764_deseq2.csv
│       ├── GSE116250_deseq2.csv
│       ├── GSE120852_deseq2.csv
│       └── GSE135055_deseq2.csv       # log2FC, p-value, padj per gene, per dataset
├── 03_scoring/
│   ├── pi_composite_scoring.xlsx      # PI-value, HSS/SSS, consistency, composite score (Excel formulas)
│   └── composite_scores.csv           # Exported flat table version, for reproducibility/version control
├── 04_intersection/
│   ├── upregulated_genes.csv          # 73 intersecting genes (70th percentile cutoff)
│   ├── downregulated_genes.csv        # 54 intersecting genes (30th percentile cutoff)
│   └── venn_diagrams/
├── 05_clustering/
│   ├── kmeans_upregulated.py          # K=4, elbow + silhouette validated
│   ├── kmeans_downregulated.py        # K=3, elbow + silhouette validated
│   └── figures/                       # elbow plots, silhouette scores, PCA cluster plots
├── 06_enrichment/
│   └── go_kegg_tables/                # g:Profiler GO/KEGG results per cluster (exported, p<0.05)
├── 07_ppi/
│   └── string_networks/               # STRING PPI networks + node degree tables (exported)
├── presentation/
└── README.md
```

## File Notes

- **Raw count matrices are not committed** — they're already public under the GEO accessions listed below; `01_datasets/data_sources.md` documents exactly which supplementary file was used per dataset.
- **Metadata files** are named `{GEO_ID}_sample_metadata.csv` and produced by `download_metadata.R`. Each dataset's `condition` column was manually standardized: stroke datasets use `"ischemic stroke"` / `"healthy"`; heart failure datasets use `"heart failure"` / `"healthy"` (see note below on GSE116250).
- **`run_deseq2.R`** is a single parameterized script (not 7 separate copies) that loops over two dataset lists — `stroke_datasets` and `heart_failure_datasets` — each pulling from one shared treatment/reference label constant, so the contrast labels stay uniform within each disease group.
- **`pi_composite_scoring.xlsx`** contains the actual formulas used (PI-value, HSS, SSS, consistency, composite score) — kept alongside a flattened `.csv` export since Excel files don't diff cleanly in git but the formulas are useful for anyone auditing the calculation.
- GO/KEGG and PPI steps were run through g:Profiler and STRING's web interfaces respectively — those folders hold exported result tables/images rather than analysis scripts.
- **Known metadata correction:** GSE116250's `condition` column originally had a labeling error (disease samples were marked `"ischemic stroke"` instead of `"heart failure"`, likely a copy-paste artifact from the stroke datasets). This was corrected before running DESeq2 — worth noting in case reviewers cross-reference the raw GEO `characteristics` field, which shows the underlying subtypes as `disease: dilated cardiomyopathy` / `disease: ischemic cardiomyopathy` / `disease: non-failing`.

## Datasets

| Stroke | Stroke n | Healthy n | &nbsp; | Heart Failure | HF n | Healthy n |
|---|---|---|---|---|---|---|
| GSE122709 | 5 | 10 | | GSE116250 | 50 | 14 |
| GSE140275 | 3 | 3 | | GSE120852 | 10 | 5 |
| GSE56267 | 7 | 6 | | GSE135055 | 7 | 6 |
| GSE178764 | 3 | 3 | | | | |

All datasets from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/).

## Pipeline

| Step | Tool | Output |
|---|---|---|
| 1. Metadata download | R (GEOquery) — `download_metadata.R` | `{GEO_ID}_sample_metadata.csv` per dataset |
| 2. Raw counts | Manual GEO download | See `data_sources.md` |
| 3. Differential expression | DESeq2 (R) — `run_deseq2.R`, merged gene-wise | log2FC, p-value per gene |
| 4. Gene scoring | Excel — PI-value (Xiao et al., 2014); HSS/SSS = avg. signed PI; Consistency = sig. count / tested count | Composite Score = 0.5×HSS + 0.5×Consistency |
| 5. Prioritization | 70th percentile (up) / 30th percentile (down) cutoff on HSS/SSS, consistency ≥ 0.5 | High-confidence up/down gene sets |
| 6. Intersection | HF ∩ IS gene sets, separately for up/down | 73 upregulated, 54 downregulated |
| 7. Clustering | Python — K-means (elbow + silhouette validated) | K=4 (up), K=3 (down) |
| 8. Enrichment | g:Profiler GO Biological Process + KEGG, p < 0.05 | Per-cluster pathway tables |
| 9. Network analysis | STRING PPI, node degree | Hub genes per direction |

## Key Results

**Upregulated hub genes (top by node degree):** AHSP, HBB, HBD, HBG1, KLF1
**Downregulated hub genes (top by node degree):** MT-ND1, MRPL3, MT-ATP6, MT-CO2, MT-ND2

| Direction | Cluster theme | Enriched pathways |
|---|---|---|
| Up | Hemoglobin / heme biosynthesis (HBB, HBD, HBG1/2, ALAS2) | Hemoglobin binding, cell development |
| Up | Vascular remodeling & inflammation (CLDN5, OSM, ITGAX, TNFRSF12A) | — |
| Down | Mitochondrial OXPHOS (ATP6, ND1, ND2, COX2) | Oxidative phosphorylation, respiratory chain complex I, ATP synthesis (p = 5.95E-08) |
| Down | Mitochondrial ribosome/translation (MRPL3, MRPL32, NDUFA12) | Mitochondrial ribosome, organellar ribosome |

## Molecular Signatures Linking HF and IS

- **Mitochondrial failure** — shared OXPHOS suppression (ATP6, ND1, COX2); acute in IS, chronic in HF
- **Shared hypoxic response** — ALAS2, KLF1, SLC4A1, hemoglobin genes upregulated in both; fetal Hb re-induction suggests convergent stress erythropoiesis
- **Common vascular remodeling** — CLDN5 (BBB), PLVAP, EGFL7, COL18A1 co-activated
- **Epigenetic/cytoskeletal suppression** — MBD2, PPL, CACYBP downregulated in both
- **Elevated inflammation** — OSM, ITGAX, TNFRSF12A, IL22RA2, SNAI2 upregulated in both
- **Actionable hub genes** — HBA2, MT-ATP6, JAK1, ALAS2, KLF1 as candidate cross-disease biomarkers/therapeutic targets

## Tools

R (GEOquery, DESeq2) · Excel (PI/composite scoring) · Python (pandas, scikit-learn, matplotlib/seaborn) · g:Profiler · STRING

## Limitations

- Small cluster sizes limit GO enrichment power
- PI-scores sensitive to outlier datasets (e.g. GSE140275); single-cohort inflation may skew rankings
- Ethnicity metadata incomplete across cohorts
- Scoring step done manually in Excel rather than a scripted pipeline — `pi_composite_scoring.xlsx` documents the formulas, but a scripted version would be more directly reproducible

## Future Directions

- Single-cell/spatial transcriptomics to resolve cell-type-specific drivers
- Mendelian randomization + prospective cohorts for causal validation of hub genes
- Convert the Excel scoring step into a reusable script for full pipeline automation

## Authors

Mehak Khemka (2509318) · Anusree P (2500146) — PGDBI 2025–26

## References

Xiao et al. (2014, *Bioinformatics*) — PI-value method · Khechaduri et al. (2013, *JACC*) · Sawicki et al. (2015, *JAHA*) · He et al. (2009, *J Cereb Blood Flow Metab*) · Chen et al. (2025, *Eur J Med Res*) · Kuswardana et al. (2025, *Jurnal Ilmiah Edutic*) — full list in presentation slides
