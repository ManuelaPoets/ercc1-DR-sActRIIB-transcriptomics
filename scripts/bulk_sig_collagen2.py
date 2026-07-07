# -*- coding: utf-8 -*-
import pandas as pd, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

GENES = ["Col3a1", "Col4a1", "Mmp2", "Timp2", "Lox"]
CONTRASTS = [("res_age_df", "age vs ctrl (aging)"),
             ("res_intv_df", "combi vs age (rescue)"),
             ("res_combi_ctrl_df", "combi vs ctrl (restoration)")]
KID = r"D:\manuela\Thesis\results\kidney_bulk\deseq2_out"
MUS = r"D:\manuela\Thesis\results\muscle_bulk\v2\deseq2_out"

# id -> symbol map from kidney table (kidney has both gene_id + gene_symbol; same Ensembl 110 annotation)
kmap = pd.read_csv(f"{KID}/res_age_df.csv")[["gene_id", "gene_symbol"]].drop_duplicates()
id2sym = dict(zip(kmap.gene_id, kmap.gene_symbol))
sym2id = {v: k for k, v in id2sym.items()}

def stars(p):
    if pd.isna(p): return "ns(NA)"
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"

def report(base, has_symbol):
    for fn, label in CONTRASTS:
        df = pd.read_csv(f"{base}/{fn}.csv")
        if not has_symbol:
            df["gene_symbol"] = df["gene_id"].map(id2sym)
        sub = df[df["gene_symbol"].isin(GENES)].set_index("gene_symbol")
        print(f"  -- {label} --")
        for g in GENES:
            if g in sub.index:
                r = sub.loc[g]
                print(f"     {g:7s} LFC={r['log2FoldChange']:+.2f}  padj={r['padj']:.2e}  {stars(r['padj'])}")
            else:
                print(f"     {g:7s} (absent)")

print("================ KIDNEY ================"); report(KID, True)
print("\n================ MUSCLE v2 ================"); report(MUS, False)
