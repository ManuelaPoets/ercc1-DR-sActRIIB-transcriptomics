#!/usr/bin/env python3
"""
Thesis supplementary analyses — PART 2
Additional receptor and MMP/TIMP markers to strengthen Discussion threads.
Run on server after Part 1.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = './thesis_supp_analysis/'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# PATHS — same as Part 1
# ============================================================
KIDNEY_ANNOTATED = "/data/bonn-epyc/projects/manuela/aging/scRNA_aging-mouse/annotated-aging-soupxed.h5ad"
MUSCLE_ANNOTATED = "/data/bonn-epyc/projects/manuela/aging/aging_muscle/annotated-muscle-soupxed2.h5ad"

KIDNEY_CONDITION_MAP = {
    'sg18': 'ctrl', 'sg20': 'age', 'sg24': 'DR', 'sg26': 'sAct', 'sg28': 'combi'
}
MUSCLE_CONDITION_MAP = {
    'rgzj1': 'ctrl', 'rgzj2': 'age', 'rgzj3': 'DR', 'rgzj4': 'sAct', 'rgzj5': 'combi'
}

# ============================================================
# LOAD DATA
# ============================================================
print("Loading kidney...")
k_adata = sc.read_h5ad(KIDNEY_ANNOTATED)
if 'sample' in k_adata.obs.columns:
    cond_col = 'sample'
elif 'orig.ident' in k_adata.obs.columns:
    cond_col = 'orig.ident'
else:
    cond_col = k_adata.obs.columns[0]
k_adata.obs['condition'] = k_adata.obs[cond_col].map(KIDNEY_CONDITION_MAP).fillna(k_adata.obs[cond_col])
k_ct_col = 'mixed_celltype' if 'mixed_celltype' in k_adata.obs.columns else 'celltype'
print(f"Kidney: {k_adata.shape}, ct col: {k_ct_col}")

print("Loading muscle...")
m_adata = sc.read_h5ad(MUSCLE_ANNOTATED)
if 'sample' in m_adata.obs.columns:
    m_cond_col = 'sample'
elif 'orig.ident' in m_adata.obs.columns:
    m_cond_col = 'orig.ident'
else:
    m_cond_col = m_adata.obs.columns[0]
m_adata.obs['condition'] = m_adata.obs[m_cond_col].map(MUSCLE_CONDITION_MAP).fillna(m_adata.obs[m_cond_col])
m_ct_col = 'mixed_celltype' if 'mixed_celltype' in m_adata.obs.columns else 'celltype'
print(f"Muscle: {m_adata.shape}, ct col: {m_ct_col}")

cond_order = ['ctrl', 'age', 'DR', 'sAct', 'combi']


# ============================================================
# ANALYSIS 6: RECEPTOR GENES ON PT (Thread 3)
# ============================================================
print("\n" + "="*60)
print("ANALYSIS 6: Receptor genes on PT-1/PT-2")
print("="*60)

receptor_genes = ['Egfr', 'Tgfbr1', 'Tgfbr2', 'Pdgfra', 'Pdgfrb',
                  'Fgfr1', 'Fgfr2', 'Met', 'Igf1r', 'Insr']

available_rec = [g for g in receptor_genes if g in k_adata.var_names]
missing_rec = [g for g in receptor_genes if g not in k_adata.var_names]
print(f"Available receptors: {available_rec}")
print(f"Missing: {missing_rec}")

if available_rec:
    # PT-1 and PT-2 split by condition
    k_pt = k_adata[k_adata.obs[k_ct_col].isin(['PT-1', 'PT-2'])].copy()
    k_pt.obs['ct_cond'] = k_pt.obs[k_ct_col].astype(str) + '_' + k_pt.obs['condition'].astype(str)
    ct_cond_order = [f'{ct}_{c}' for ct in ['PT-1', 'PT-2'] for c in cond_order]
    ct_cond_order = [x for x in ct_cond_order if x in k_pt.obs['ct_cond'].unique()]

    sc.pl.dotplot(k_pt, var_names=available_rec, groupby='ct_cond',
                  categories_order=ct_cond_order, standard_scale='var',
                  show=False, figsize=(12, 6))
    plt.title('PT-1 vs PT-2: receptor expression across conditions')
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}07_PT_receptors_by_condition.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: 07_PT_receptors_by_condition.png")

    # Also show across ALL kidney cell types for context
    sc.pl.dotplot(k_adata, var_names=available_rec, groupby=k_ct_col,
                  standard_scale='var', show=False, figsize=(14, 8))
    plt.title('Receptor genes across all kidney cell types')
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}08_receptors_all_kidney_celltypes.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: 08_receptors_all_kidney_celltypes.png")

    # Export CSV
    rec_expr = pd.DataFrame()
    for ct in ['PT-1', 'PT-2', 'FIB', 'vSMC/MC', 'PODO', 'IMM']:
        for cond in cond_order:
            mask = (k_adata.obs[k_ct_col] == ct) & (k_adata.obs['condition'] == cond)
            subset = k_adata[mask]
            if subset.shape[0] > 0:
                avail_here = [g for g in available_rec if g in subset.var_names]
                if hasattr(subset.X, 'toarray'):
                    means = pd.DataFrame(subset[:, avail_here].X.toarray(),
                                         columns=avail_here).mean()
                else:
                    means = pd.DataFrame(subset[:, avail_here].X,
                                         columns=avail_here).mean()
                means['celltype'] = ct
                means['condition'] = cond
                means['n_cells'] = subset.shape[0]
                rec_expr = pd.concat([rec_expr, means.to_frame().T], ignore_index=True)

    rec_expr.to_csv(f'{OUTPUT_DIR}07_receptor_expression.csv', index=False)
    print("Saved: 07_receptor_expression.csv")


# ============================================================
# ANALYSIS 7: MMP/TIMP MARKERS (Thread 1 — ECM remodelling)
# ============================================================
print("\n" + "="*60)
print("ANALYSIS 7: MMP/TIMP ECM remodelling markers")
print("="*60)

mmp_genes = ['Mmp2', 'Mmp3', 'Mmp9', 'Mmp14',
             'Timp1', 'Timp2', 'Timp3',
             'Ctgf', 'Serpine1', 'Lox', 'Loxl2']

available_mmp = [g for g in mmp_genes if g in k_adata.var_names]
missing_mmp = [g for g in mmp_genes if g not in k_adata.var_names]
print(f"Available MMP/TIMP: {available_mmp}")
print(f"Missing: {missing_mmp}")

# Try alternate name for Ctgf
if 'Ctgf' in missing_mmp and 'Ccn2' in k_adata.var_names:
    available_mmp.append('Ccn2')
    print("  Found alternate: Ccn2 for Ctgf")

if available_mmp:
    # FIB + vSMC/MC + PT-1 + PT-2 + IMM across conditions
    ecm_cts = ['FIB', 'vSMC/MC', 'PT-1', 'PT-2', 'IMM', 'PODO']
    ecm_cts_avail = [ct for ct in ecm_cts if ct in k_adata.obs[k_ct_col].unique()]
    k_ecm = k_adata[k_adata.obs[k_ct_col].isin(ecm_cts_avail)].copy()
    k_ecm.obs['ct_cond'] = k_ecm.obs[k_ct_col].astype(str) + '_' + k_ecm.obs['condition'].astype(str)

    ct_cond_order = [f'{ct}_{c}' for ct in ecm_cts_avail for c in cond_order]
    ct_cond_order = [x for x in ct_cond_order if x in k_ecm.obs['ct_cond'].unique()]

    sc.pl.dotplot(k_ecm, var_names=available_mmp, groupby='ct_cond',
                  categories_order=ct_cond_order, standard_scale='var',
                  show=False, figsize=(12, 12))
    plt.title('Kidney: MMP/TIMP ECM remodelling markers across cell types and conditions')
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}09_kidney_MMP_TIMP_markers.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: 09_kidney_MMP_TIMP_markers.png")

    # Export MMP/TIMP mean expression for FIB across conditions
    mmp_expr = pd.DataFrame()
    for ct in ecm_cts_avail:
        for cond in cond_order:
            mask = (k_adata.obs[k_ct_col] == ct) & (k_adata.obs['condition'] == cond)
            subset = k_adata[mask]
            if subset.shape[0] > 0:
                avail_here = [g for g in available_mmp if g in subset.var_names]
                if hasattr(subset.X, 'toarray'):
                    means = pd.DataFrame(subset[:, avail_here].X.toarray(),
                                         columns=avail_here).mean()
                else:
                    means = pd.DataFrame(subset[:, avail_here].X,
                                         columns=avail_here).mean()
                means['celltype'] = ct
                means['condition'] = cond
                means['n_cells'] = subset.shape[0]
                mmp_expr = pd.concat([mmp_expr, means.to_frame().T], ignore_index=True)

    mmp_expr.to_csv(f'{OUTPUT_DIR}09_MMP_TIMP_expression.csv', index=False)
    print("Saved: 09_MMP_TIMP_expression.csv")

    # MMP2/TIMP2 ratio in FIB (enzymatic remodelling capacity)
    if 'Mmp2' in available_mmp and 'Timp2' in available_mmp:
        ratio_data = []
        for cond in cond_order:
            mask = (k_adata.obs[k_ct_col] == 'FIB') & (k_adata.obs['condition'] == cond)
            subset = k_adata[mask]
            if subset.shape[0] > 0:
                if hasattr(subset.X, 'toarray'):
                    mmp2 = subset[:, 'Mmp2'].X.toarray().mean()
                    timp2 = subset[:, 'Timp2'].X.toarray().mean()
                else:
                    mmp2 = subset[:, 'Mmp2'].X.mean()
                    timp2 = subset[:, 'Timp2'].X.mean()
                ratio = mmp2 / timp2 if timp2 > 0 else np.nan
                ratio_data.append({'condition': cond, 'Mmp2_mean': mmp2,
                                   'Timp2_mean': timp2, 'ratio_MMP2_TIMP2': ratio,
                                   'n_cells': subset.shape[0]})
        ratio_df = pd.DataFrame(ratio_data)
        ratio_df.to_csv(f'{OUTPUT_DIR}10_kidney_FIB_MMP2_TIMP2_ratio.csv', index=False)
        print("Saved: 10_kidney_FIB_MMP2_TIMP2_ratio.csv")
        print(ratio_df)


# ============================================================
# ANALYSIS 8: PDGF pathway genes (Thread 1 — colleague suggestion)
# ============================================================
print("\n" + "="*60)
print("ANALYSIS 8: PDGF pathway genes")
print("="*60)

pdgf_genes = ['Pdgfa', 'Pdgfb', 'Pdgfc', 'Pdgfd',
              'Pdgfra', 'Pdgfrb']

available_pdgf = [g for g in pdgf_genes if g in k_adata.var_names]
print(f"Available PDGF genes: {available_pdgf}")

if available_pdgf:
    # FIB + vSMC/MC + EC subtypes + PODO
    pdgf_cts = ['FIB', 'vSMC/MC', 'PODO', 'EC-1(gEC)', 'EC-2', 'EC-3', 'PT-1', 'PT-2', 'IMM']
    pdgf_cts_avail = [ct for ct in pdgf_cts if ct in k_adata.obs[k_ct_col].unique()]
    k_pdgf = k_adata[k_adata.obs[k_ct_col].isin(pdgf_cts_avail)].copy()
    k_pdgf.obs['ct_cond'] = k_pdgf.obs[k_ct_col].astype(str) + '_' + k_pdgf.obs['condition'].astype(str)

    ct_cond_order = [f'{ct}_{c}' for ct in pdgf_cts_avail for c in cond_order]
    ct_cond_order = [x for x in ct_cond_order if x in k_pdgf.obs['ct_cond'].unique()]

    sc.pl.dotplot(k_pdgf, var_names=available_pdgf, groupby='ct_cond',
                  categories_order=ct_cond_order, standard_scale='var',
                  show=False, figsize=(10, 14))
    plt.title('Kidney: PDGF ligands and receptors across cell types and conditions')
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}11_kidney_PDGF_pathway.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: 11_kidney_PDGF_pathway.png")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*60)
print("PART 2 — ALL NEW OUTPUTS:")
print("="*60)
new_files = [f for f in sorted(os.listdir(OUTPUT_DIR)) if f.startswith(('07', '08', '09', '10', '11'))]
for f in new_files:
    size = os.path.getsize(f'{OUTPUT_DIR}{f}')
    print(f"  {f} ({size/1024:.1f} KB)")

print(f"\nAll files saved to: {OUTPUT_DIR}")
print("Upload CSVs and PNGs to Claude for review.")
