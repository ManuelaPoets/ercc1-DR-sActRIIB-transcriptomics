# -*- coding: utf-8 -*-
import pandas as pd, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

GENES = ["Col3a1", "Col4a1", "Mmp2", "Timp2", "Lox"]
CONTRASTS = [("res_age_df", "age vs ctrl (aging)"),
             ("res_intv_df", "combi vs age (rescue)"),
             ("res_combi_ctrl_df", "combi vs ctrl (restoration)")]
BASES = {"KIDNEY": r"D:\manuela\Thesis\results\kidney_bulk\deseq2_out",
         "MUSCLE v2": r"D:\manuela\Thesis\results\muscle_bulk\v2\deseq2_out"}

def stars(p):
    if pd.isna(p): return "ns(NA)"
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

for tissue, base in BASES.items():
    print(f"\n================ {tissue} ================")
    for fn, label in CONTRASTS:
        df = pd.read_csv(f"{base}/{fn}.csv")
        sym = "gene_symbol" if "gene_symbol" in df.columns else df.columns[1]
        sub = df[df[sym].isin(GENES)].set_index(sym)
        print(f"  -- {label} --")
        for g in GENES:
            if g in sub.index:
                r = sub.loc[g]
                print(f"     {g:7s} LFC={r['log2FoldChange']:+.2f}  padj={r['padj']:.2e}  {stars(r['padj'])}")
            else:
                print(f"     {g:7s} (not in table / filtered out)")
