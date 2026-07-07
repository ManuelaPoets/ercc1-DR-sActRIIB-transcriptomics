# Supplementary Tables (S2–S11)

Curated data tables supporting the dissertation
*"Transcriptomic Rescue of Accelerated Aging by Combined Dietary Restriction and
Activin Signaling Inhibition in Ercc1Δ/− Kidney and Muscle"* (M. W. Poets).

Table S1 (sample metadata) is printed in the thesis appendix. Raw and processed
sequencing data are on NCBI GEO under accession **GSE268971**.

## Provenance (important)

All **kidney** rescue tables here are the authoritative **post-shrinkage** export
(ashr lfcShrink applied to the combi-vs-ctrl contrast; rescued-to-normal = 2,422).
Earlier pre-shrinkage kidney tables (rescued-to-normal = 1,961) are **stale and are
not included**. Muscle (v2, STAR+Salmon/GRCm39) tables are current as exported.

## Folders

| Folder | Contents |
|--------|----------|
| `kidney_shrunk_authoritative/` | All kidney rescue tables (post-shrink): contrasts, rescued-to-normal (2,422), significant rescued-to-normal (277), non-rescued (1,169), combined reversal (4,645), intervention-impact / driver categories (2,422; DR 652 / sAct 156 / combi-only 109 / dual 1,505), per-cell-type markers (1,142 hits / 423 genes; PT-2 132, FIB 115, vSMC/MC 92, EC 225), per-cell-type GO:BP, bulk enrichment (RTN GO:BP 765, DR 257), and kidney CellChat L–R tables |
| `muscle_v2/` | Muscle rescue tables (rescued-to-normal 491; drivers 211/20/14/246), enrichment, effect-comparison, CellChat L–R |
| `cross_tissue/` | 95-gene kidney∩muscle rescued core (`Common_rescued_genes_kidney_muscle.csv`; r = 0.658, 81/95 same-direction), plus muscle directional GO:BP splits |
| `enrichment_DAVID/` | DAVID functional-annotation clusters, kidney (post-shrink) and muscle |

## Index (S-numbers → files)

| ID | Content | Location |
|----|---------|----------|
| S2 | GO:BP enrichment (bulk) | `kidney_shrunk_authoritative/enrichment/`, `muscle_v2/2_Enrichment_bulk&sn/`, `cross_tissue/Muscle_*` |
| S3 | DAVID clusters | `enrichment_DAVID/` |
| S4 | Non-rescued gene sets | `kidney_shrunk_authoritative/kidney_non_rescued_aged_DEGs.csv`, `muscle_v2/1_DEGs_tables/*_non_rescued_aged_DEGs.csv` |
| S5 | Full aging DEG results | `kidney_shrunk_authoritative/kidney_res_age_results_with_symbols.csv`, `muscle_v2/1_DEGs_tables/*_res_age_*` |
| S6 | Rescued-to-normal lists | `kidney_shrunk_authoritative/kidney_rescued_to_normal.csv`, `muscle_v2/1_DEGs_tables/*_rescued_to_normal.csv` |
| S7 | Intervention impact / drivers | `*/intervention_impact/` and `muscle_v2/3_Effect_comparison_tables/` |
| S8 | Per-cell-type reversal markers + GO:BP | `kidney_shrunk_authoritative/celltype_markers/`, `celltype_pathways/`, `muscle_v2/…` |
| S9 | CellChat L–R interaction tables | `kidney_shrunk_authoritative/cellchat_outputs_kidney/`, `muscle_v2/5_Cellchat_LR_analysis/` |
| S10 | 95 common rescued genes | `cross_tissue/Common_rescued_genes_kidney_muscle.csv` |
| S11 | Exploratory snRNA-seq cell-type DEG tables | `kidney_shrunk_authoritative/`, `muscle_v2/1_DEGs_tables/` |
