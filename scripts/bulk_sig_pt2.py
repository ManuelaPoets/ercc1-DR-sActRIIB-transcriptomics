# -*- coding: utf-8 -*-
import pandas as pd, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
GENES = ["Havcr1", "Cdkn1a", "Slc34a1", "Vcam1", "Aqp1"]
CONTRASTS = [("res_age_df", "age vs ctrl (aging)"),
             ("res_intv_df", "combi vs age (rescue)"),
             ("res_combi_ctrl_df", "combi vs ctrl (restoration)")]
KID = r"D:\manuela\Thesis\results\kidney_bulk\deseq2_out"
def stars(p):
    if pd.isna(p): return "ns(NA)"
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"
print("================ KIDNEY (PT-2 markers, whole-tissue bulk) ================")
for fn, label in CONTRASTS:
    df = pd.read_csv(f"{KID}/{fn}.csv")
    sub = df[df["gene_symbol"].isin(GENES)].set_index("gene_symbol")
    print(f"  -- {label} --")
    for g in GENES:
        if g in sub.index:
            r = sub.loc[g]
            print(f"     {g:8s} LFC={r['log2FoldChange']:+.2f}  padj={r['padj']:.2e}  {stars(r['padj'])}")
        else:
            print(f"     {g:8s} (absent)")
