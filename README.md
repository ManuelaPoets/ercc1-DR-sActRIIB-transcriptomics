# Transcriptomic rescue of accelerated aging by dietary restriction and activin signaling inhibition in Ercc1Δ/− kidney and muscle

Analysis code accompanying the doctoral thesis of **Manuela Weronika Poets**
(Institute for Medical Systems Bioinformatics / III. Medizinische Klinik und
Poliklinik, Universitätsklinikum Hamburg-Eppendorf; supervisors Prof. Dr. med.
Tobias B. Huber and Prof. Dr. Stefan Bonn).

This repository contains the custom R and Python analysis code used to process
bulk RNA-seq and single-nucleus RNA-seq (snRNA-seq) data from the Ercc1Δ/−
progeroid mouse model and to run the direction-aware rescue framework, pathway
enrichment, cell-type mapping, cell–cell communication (CellChat), and
cross-tissue integration reported in the thesis.

> **This repository holds code only.** No raw or processed sequencing data are
> included (see *Data availability*). Thesis-document tooling (figure embedding,
> comment/track-change scripts) is intentionally excluded.

## Repository structure

```
notebooks/            Core analysis notebooks (the primary pipeline)
  DESeq2_bulk_sn/      Bulk DESeq2 + direction-aware rescue framework (kidney, muscle v2)
  sn_analysis/         snRNA-seq: SoupX, integration, annotation, DE, module scores
  aging_interorgan/    Cross-tissue merge + inter-organ CellChat
  Discussion_analysis/ PT1/PT2 re-annotation, collagen/laminin, supplementary analyses
  Google_colab/        Results-figure plotting
scripts/              R/Python analysis + verification scripts and figure generators
  (DESeq2 reproduction, rescue counts, CellChat inspection, per-figure discussion
   analyses, and val_/verify_ data-verification checks)
envs/                 Pinned software versions (Python, R, alignment tools)
```

A list of files deliberately excluded from the original working folder (document
tooling and superseded iterations) is in `_MANIFEST_excluded.txt`.

## Pipeline overview (maps to thesis Methods §2)

1. **Bulk RNA-seq** — nf-core/rnaseq v3.9 (STAR 2.7.10a + Salmon 1.5.2, GRCm39/
   Ensembl 110); import with tximport; DESeq2 1.46.0 with ashr log2FC shrinkage
   (`notebooks/DESeq2_bulk_sn/`, `scripts/*deseq2*.R`).
2. **Direction-aware rescue framework** — aging signature → rescued-to-normal
   (Stage 1: sign reversal + |log2FC combi vs ctrl| < 0.5) → significant subset
   (Stage 2) → DR/sActRIIB driver classification and restoration scores
   (`notebooks/DESeq2_bulk_sn/`, `scripts/reproduce_rescue_exact.R`,
   `scripts/muscle_v2_rescue_counts.R`).
3. **snRNA-seq** — SoupX 1.6.2 → Scanpy 1.9.6 (Scrublet, QC, Harmony, Leiden,
   UMAP) → expert cell-type annotation (`notebooks/sn_analysis/`).
4. **Bulk–single-cell integration** — Seurat 5.4.0 AddModuleScore projection of
   bulk-derived gene sets; reversal-gene-marker overlap; per-cell-type enrichGO.
5. **Cell–cell communication** — CellChat 1.6.1, intra-tissue and inter-organ
   (`notebooks/aging_interorgan/`, `scripts/*cellchat*.R`).
6. **Cross-tissue integration** — shared rescued-gene core, directional
   consensus, Enrichr/GO enrichment (`notebooks/aging_interorgan/3_cross-tissue_analysis.ipynb`).

## Environments

See `envs/`. Core pins: Python 3.10 (scanpy 1.9.6, anndata 0.7.8, numpy<2,
pandas<2); R 4.4.3 (DESeq2 1.46.0, Seurat 5.4.0, CellChat 1.6.1). Known
workaround: call `harmonypy::run_harmony()` directly rather than
`sc.external.pp.harmony_integrate()`.

## Data availability

- **Bulk RNA-seq** (kidney + muscle): NCBI GEO **GSE268971** (Vermeij et al., 2024),
  re-processed here from raw reads.
- **snRNA-seq** (kidney + muscle): to be deposited upon publication.
- Reference genome: *Mus musculus* GRCm39, Ensembl release 110.

## Citation

If you use this code, please cite the thesis (Poets, M. W., Universität Hamburg,
2026) and the primary dataset (Vermeij et al., *J Cachexia Sarcopenia Muscle*
2024; GSE268971).

## License

MIT — see `LICENSE`.
