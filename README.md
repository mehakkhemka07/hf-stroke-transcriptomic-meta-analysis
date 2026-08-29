# Consensus Signature Markers Interconnecting Heart Failure and Ischemic Stroke

A multi-dataset transcriptomic meta-analysis identifying shared gene expression signatures between heart failure (HF) and ischemic stroke (IS) across 7 GEO datasets, using DESeq2, a Priority Index (PI) scoring system, K-means clustering, GO/KEGG enrichment, and STRING PPI network analysis to prioritize cross-disease hub genes.

## Repository Structure

```
hf-stroke-transcriptomic-meta-analysis/
├── 01_datasets/
│   ├── stroke/                        # GSE122709, GSE140275, GSE56267, GSE178764
│   └── heart_failure/                 # GSE116250, GSE120852, GSE135055
├── 02_deseq2/
│   ├── scripts/                       # Per-dataset DESeq2 pipelines (R)
│   └── results/                       # log2FC, p-value, padj per gene per dataset
├── 03_scoring/
│   ├── pi_score_calculation.py        # PI = sign(log2FC) × |log2FC| × -log10(p-value)
│   └── composite_scores.csv           # HSS/SSS, consistency, composite score per gene
├── 04_intersection/
│   ├── upregulated_genes.csv          # 73 intersecting genes
│   ├── downregulated_genes.csv        # 54 intersecting genes
│   └── venn_diagrams/
├── 05_clustering/
│   ├── kmeans_upregulated.py          # K=4
│   ├── kmeans_downregulated.py        # K=3
│   └── figures/                       # elbow plots, silhouette scores, PCA clusters
├── 06_enrichment/
│   └── go_kegg_tables/                # g:Profiler results per cluster
├── 07_ppi/
│   └── string_networks/               # PPI networks + node degree tables
├── presentation/
└── README.md
```

## Datasets

| Stroke | Stroke n | Healthy n | &nbsp; | Heart Failure | HF n | Healthy n |
|---|---|---|---|---|---|---|
| GSE122709 | 5 | 10 | | GSE116250 | 50 | 14 |
| GSE140275 | 3 | 3 | | GSE120852 | 10 | 5 |
| GSE56267 | 7 | 6 | | GSE135055 | 7 | 6 |
| GSE178764 | 3 | 3 | | | | |

All datasets from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/).

## Pipeline

| Step | Method | Output |
|---|---|---|
| 1. Differential expression | DESeq2 per dataset, merged gene-wise | log2FC, p-value per gene |
| 2. Gene scoring | PI-value (Xiao et al., 2014); HSS/SSS = avg. signed PI; Consistency = sig. count / tested count | Composite Score = 0.5×HSS + 0.5×Consistency |
| 3. Prioritization | 70th percentile (up) / 30th percentile (down) cutoff on HSS/SSS, consistency ≥ 0.5 | High-confidence up/down gene sets |
| 4. Intersection | HF ∩ IS gene sets, separately for up/down | 73 upregulated, 54 downregulated |
| 5. Clustering | K-means (elbow + silhouette validated) | K=4 (up), K=3 (down) |
| 6. Enrichment | g:Profiler GO Biological Process + KEGG, p < 0.05 | Per-cluster pathway tables |
| 7. Network analysis | STRING PPI, node degree | Hub genes per direction |

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

R (DESeq2) · Python (pandas, scikit-learn, matplotlib/seaborn) · g:Profiler · STRING · Excel

## Limitations

- Small cluster sizes limit GO enrichment power
- PI-scores sensitive to outlier datasets (e.g. GSE140275); single-cohort inflation may skew rankings
- Ethnicity metadata incomplete across cohorts

## Future Directions

- Single-cell/spatial transcriptomics to resolve cell-type-specific drivers
- Mendelian randomization + prospective cohorts for causal validation of hub genes

## Authors

Mehak Khemka (2509318) · Anusree P (2500146) — PGDBI 2025–26

## References

Xiao et al. (2014, *Bioinformatics*) — PI-value method · Khechaduri et al. (2013, *JACC*) · Sawicki et al. (2015, *JAHA*) · He et al. (2009, *J Cereb Blood Flow Metab*) · Chen et al. (2025, *Eur J Med Res*) · Kuswardana et al. (2025, *Jurnal Ilmiah Edutic*) — full list in presentation slides
