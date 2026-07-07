import h5py, numpy as np
from scipy.sparse import csr_matrix

path = r"D:\manuela\Thesis\results\sn\annotated-aging-soupxed.h5ad"
f = h5py.File(path, "r")
g = f["raw"]["X"]; shape = tuple(g.attrs["shape"])
X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).tocsc()
genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["raw"]["var"]["_index"][:]])
gidx = {n: int(np.where(genes==n)[0][0]) for n in ["Col3a1","Col4a1","Mmp2","Timp2"]}
print("raw/X max:", round(float(X.max()),2), "(log-norm if <15)")
def cat(col):
    o=f["obs"][col]; cats=np.array([c.decode() if isinstance(c,bytes) else c for c in o["categories"][:]]); return cats[o["codes"][:]]
celltype = cat("celltype"); condition = cat("sample")
f.close()

conds = ["ctrl","age","DR","sAct","combi"]
cols = {gn: np.asarray(X[:, gidx[gn]].todense()).ravel() for gn in gidx}
ct="FIB"
print("FIB n per condition:", {c:int(((celltype==ct)&(condition==c)).sum()) for c in conds})
def stat(gn,c):
    m=(celltype==ct)&(condition==c); v=cols[gn][m]
    return dict(meanLog=v.mean(), det=100*np.mean(v>0), logCP=np.log1p(np.expm1(v).mean()))
for num,den in [("Col3a1","Col4a1"),("Mmp2","Timp2")]:
    print(f"\n== {num}/{den} (kidney FIB) ==")
    for c in conds:
        sn=stat(num,c); sd=stat(den,c)
        print(f"  {c:5s}: {num}(meanLog={sn['meanLog']:.3f} det={sn['det']:.0f}%) {den}(meanLog={sd['meanLog']:.3f} det={sd['det']:.0f}%) "
              f"| ratio_meanLog={sn['meanLog']/sd['meanLog']:.2f}  ratio_logCP={sn['logCP']/sd['logCP']:.2f}")
print("""
THESIS (4.3) kidney FIB: Col3a1/Col4a1 ctrl2.90 age3.31 DR3.15 sAct2.89 combi1.16
  (Col4a1 0.80->0.73 ; Col3a1 2.31->2.41 ; combi Col3a1 0.87 Col4a1 0.75 ; DR Col3a1 1.60)
  Mmp2/Timp2 ctrl0.23 age1.41 combi0.83 (Mmp2 0.272->1.213 ; Timp2 1.17->0.86 ; combi Mmp2 0.52 Timp2 0.62)
""")
