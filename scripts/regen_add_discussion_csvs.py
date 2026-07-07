#!/usr/bin/env python3
"""Regenerate the 7 mismatched CSVs in results/add_discussion_analysis/.

Bug in originals: ratio/expression means were taken from the SCALED X (z-scored
HVG matrix) instead of the lognorm values, producing ratios 4-10x off. Detection
files likewise need lognorm-based detection. This recomputes everything from
raw.X (full-gene lognorm). Detection % = 100*mean(raw.X>0), identical to
100*mean(counts>0) because log1p(CPM) is zero iff the count is zero.

Outputs go to results/add_discussion_analysis/corrected/ (originals untouched).
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

ROOT = "/media/mpoets/Expansion/manuela/Thesis"
OUT = os.path.join(ROOT, "results/add_discussion_analysis/corrected")
os.makedirs(OUT, exist_ok=True)

COND_ORDER = ["ctrl", "age", "DR", "sAct", "combi"]
MUSCLE_MAP = {"rgzj1": "ctrl", "rgzj2": "age", "rgzj3": "sAct",
              "rgzj4": "DR", "rgzj5": "combi"}


def get_raw(adata):
    """Return (raw_X csr, var_names list)."""
    X = adata.raw.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    return sp.csr_matrix(X), list(adata.raw.var_names)


def col(rawX, vn_index, gene):
    j = vn_index[gene]
    return np.asarray(rawX[:, j].todense()).ravel()


def ratio_table(adata, rawX, vn_index, celltype, g1, g2, extra_cols=None):
    ct = adata.obs["celltype"].astype(str).values
    cond = adata.obs["__cond"].values
    rows = []
    for c in COND_ORDER:
        mask = (ct == celltype) & (cond == c)
        n = int(mask.sum())
        m1 = float(col(rawX, vn_index, g1)[mask].mean()) if n else np.nan
        m2 = float(col(rawX, vn_index, g2)[mask].mean()) if n else np.nan
        ratio = m1 / m2 if m2 else np.nan
        row = {"condition": c, f"{g1}_mean": m1, f"{g2}_mean": m2,
               "ratio": ratio, "n_cells": n}
        if extra_cols:
            row.update(extra_cols)
        rows.append(row)
    return pd.DataFrame(rows)


def detection_table(adata, rawX, vn_index, celltype, markers):
    ct = adata.obs["celltype"].astype(str).values
    cond = adata.obs["__cond"].values
    rows = []
    for c in COND_ORDER:
        mask = (ct == celltype) & (cond == c)
        n = int(mask.sum())
        row = {"condition": c, "n_cells": n}
        for m in markers:
            if m not in vn_index:
                row[f"{m}_pct"] = np.nan
                row[f"{m}_mean"] = np.nan
                continue
            v = col(rawX, vn_index, m)[mask]
            row[f"{m}_pct"] = float(100.0 * (v > 0).mean()) if n else np.nan
            row[f"{m}_mean"] = float(v.mean()) if n else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------- KIDNEY ----------------
print("loading kidney...")
k = ad.read_h5ad(os.path.join(ROOT, "original_data/sn/annotated-aging-soupxed.h5ad"))
k.obs["__cond"] = k.obs["sample"].astype(str).values
kraw, kvn = get_raw(k)
kidx = {g: i for i, g in enumerate(kvn)}

t03 = ratio_table(k, kraw, kidx, "FIB", "Col3a1", "Col4a1")
t10 = ratio_table(k, kraw, kidx, "FIB", "Mmp2", "Timp2")
IMM_MARKERS = ["Adgre1", "Cd14", "Cd68", "Cd80", "Clec9a", "Csf1r", "Fcgr3",
               "H2-Aa", "H2-Ab1", "H2-Eb1", "Kit", "Ms4a1", "Ptprc"]
t23 = detection_table(k, kraw, kidx, "IMM", IMM_MARKERS)

# validation: reproduce canonical-22 PT-1/PT-2 detection for a couple markers
val_pt = detection_table(k, kraw, kidx, "PT-1", ["Vcam1", "Lrp2", "Havcr1"])

# ---------------- MUSCLE ----------------
print("loading muscle...")
m = ad.read_h5ad(os.path.join(ROOT, "original_data/sn/annotated-muscle-soupxed.h5ad"))
m.obs["__cond"] = m.obs["sample"].astype(str).map(MUSCLE_MAP).values
mraw, mvn = get_raw(m)
midx = {g: i for i, g in enumerate(mvn)}

t13 = ratio_table(m, mraw, midx, "FAPs", "Col3a1", "Col4a1")
t15 = ratio_table(m, mraw, midx, "FAPs", "Mmp2", "Timp2")
MAC_MARKERS = ["Adgre1", "Cd68", "Cd80", "Clec9a", "Csf1r", "Fcgr3",
               "H2-Aa", "H2-Ab1", "H2-Eb1", "Kit", "Ptprc"]
t25 = detection_table(m, mraw, midx, "Macrophages", MAC_MARKERS)

# 19 cross-tissue collagen = kidney FIB + muscle FAPs stacked
a19 = t03.copy(); a19["tissue"] = "Kidney (FIB)"; a19["celltype"] = "FIB"
b19 = t13.copy(); b19["tissue"] = "Muscle (FAPs)"; b19["celltype"] = "FAPs"
t19 = pd.concat([a19, b19], ignore_index=True)

# add celltype col to muscle ratio files (originals had it)
t13["celltype"] = "FAPs"; t15["celltype"] = "FAPs"

# ---------------- WRITE ----------------
writes = {
    "03_kidney_FIB_collagen_ratio.csv": t03,
    "10_kidney_FIB_MMP2_TIMP2_ratio.csv": t10,
    "13_muscle_FAP_collagen_ratio.csv": t13,
    "15_muscle_FAP_MMP2_TIMP2_ratio.csv": t15,
    "19_cross_tissue_collagen_ratio.csv": t19,
    "23_kidney_IMM_comprehensive_detection.csv": t23,
    "25_muscle_Macrophages_detection.csv": t25,
}
for fn, df in writes.items():
    df.to_csv(os.path.join(OUT, fn), index=False)

# ---------------- VALIDATION ----------------
print("\n=== VALIDATION vs handoff expected (lognorm ratios) ===")
def rget(df, c): return float(df.loc[df.condition == c, "ratio"].iloc[0])
print(f"03 FIB Col3a1/Col4a1  ctrl={rget(t03,'ctrl'):.3f} (exp 0.686) "
      f"sAct={rget(t03,'sAct'):.3f} (exp 0.627) combi={rget(t03,'combi'):.3f} (exp 0.322)")
print(f"10 FIB Mmp2/Timp2     ctrl={rget(t10,'ctrl'):.3f} (exp 0.087) "
      f"combi={rget(t10,'combi'):.3f} (exp 0.230)")
print(f"15 FAP Mmp2/Timp2     ctrl={rget(t15,'ctrl'):.3f} (exp 0.843) "
      f"combi={rget(t15,'combi'):.3f} (exp 1.226)")
print(f"13 FAP Col3a1/Col4a1  ctrl={rget(t13,'ctrl'):.3f} (no exp given)")
print("\n=== canonical-22 reproduction check (PT-1 detection %) ===")
print(val_pt[["condition", "n_cells", "Vcam1_pct", "Lrp2_pct", "Havcr1_pct"]].to_string(index=False))
print("  (22 canonical PT-1: Vcam1_pct ctrl=0.78 age=8.15; Lrp2_pct ctrl=99.06; Havcr1_pct ctrl=1.70)")
print("\n=== 23 IMM MHC-II detection (thesis-cited) ===")
print(t23[["condition", "n_cells", "H2-Aa_pct", "H2-Ab1_pct", "H2-Eb1_pct", "Adgre1_pct", "Csf1r_pct"]].to_string(index=False))
print("\n=== 25 Macrophage MHC-II detection (thesis-cited) ===")
print(t25[["condition", "n_cells", "H2-Aa_pct", "H2-Ab1_pct", "H2-Eb1_pct", "Adgre1_pct", "Csf1r_pct"]].to_string(index=False))
print("\nwrote", len(writes), "CSVs to", OUT)
