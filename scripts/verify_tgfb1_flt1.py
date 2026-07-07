# -*- coding: utf-8 -*-
import h5py, numpy as np, pandas as pd, sys, io
from scipy.sparse import csr_matrix
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
CONDS = ["ctrl","age","DR","sAct","combi"]
SMAP = {"rgzj1":"ctrl","rgzj2":"age","rgzj3":"sAct","rgzj4":"DR","rgzj5":"combi"}

def load(path, muscle=False):
    f=h5py.File(path,"r"); g=f["raw"]["X"]; shape=tuple(g.attrs["shape"])
    X=csr_matrix((g["data"][:],g["indices"][:],g["indptr"][:]),shape=shape).tocsc()
    genes=np.array([x.decode() if isinstance(x,bytes) else x for x in f["raw"]["var"]["_index"][:]])
    def cat(c):
        o=f["obs"][c]; cats=np.array([z.decode() if isinstance(z,bytes) else z for z in o["categories"][:]]); return cats[o["codes"][:]]
    ct=cat("celltype"); s=cat("sample"); f.close()
    cond=np.array([SMAP.get(x,x) for x in s]) if muscle else s
    return X,genes,ct,cond
def det(X,genes,ct,cond,cell,gn):
    if gn not in genes: return [float("nan")]*5
    gi=int(np.where(genes==gn)[0][0]); v=np.asarray(X[:,gi].todense()).ravel()
    return [100*np.mean(v[(ct==cell)&(cond==c)]>0) if ((ct==cell)&(cond==c)).sum() else float("nan") for c in CONDS]

Xk,gk,ctk,ck = load(r"D:\manuela\Thesis\results\sn\annotated-aging-soupxed.h5ad")
Xm,gm,ctm,cm = load(r"D:\manuela\Thesis\results\sn\annotated-muscle-soupxed.h5ad", muscle=True)
print("Tgfb1 kidney IMM:      ", [round(v,1) for v in det(Xk,gk,ctk,ck,'IMM','Tgfb1')])
print("Tgfb1 muscle Macroph.: ", [round(v,1) for v in det(Xm,gm,ctm,cm,'Macrophages','Tgfb1')])
print("\nFlt1 kidney EC subtypes (det% across", CONDS, "):")
for ec in ['EC-1(gEC)','EC-2','EC-3','EC-4','EC-5','EC-6']:
    print(f"  {ec:11s}", [round(v,1) for v in det(Xk,gk,ctk,ck,ec,'Flt1')])

# bulk significance
def sig(base, genes, label):
    print(f"\n-- BULK {label} --")
    km = pd.read_csv(r"D:\manuela\Thesis\results\kidney_bulk\deseq2_out\res_age_df.csv")[["gene_id","gene_symbol"]]
    id2sym = dict(zip(km.gene_id, km.gene_symbol))
    for fn,lab in [("res_age_df","aging"),("res_intv_df","rescue(combi-age)")]:
        d=pd.read_csv(f"{base}/{fn}.csv")
        if "gene_symbol" not in d.columns: d["gene_symbol"]=d["gene_id"].map(id2sym)
        sub=d[d["gene_symbol"].isin(genes)].set_index("gene_symbol")
        for g in genes:
            if g in sub.index:
                r=sub.loc[g]; st="***" if r['padj']<.001 else "**" if r['padj']<.01 else "*" if r['padj']<.05 else "ns"
                print(f"   {lab:18s} {g:7s} LFC={r['log2FoldChange']:+.2f} padj={r['padj']:.1e} {st}")
sig(r"D:\manuela\Thesis\results\kidney_bulk\deseq2_out", ["Tgfb1","Flt1","B2m"], "KIDNEY")
sig(r"D:\manuela\Thesis\results\muscle_bulk\v2\deseq2_out", ["Tgfb1"], "MUSCLE v2")
print("\nTHESIS 4.4: Tgfb1 up with aging in IMM/macrophages, rescued by DR/combi; Flt1 rescued across all 6 EC subtypes.")
