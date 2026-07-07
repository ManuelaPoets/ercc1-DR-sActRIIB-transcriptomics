import h5py, numpy as np
from scipy.sparse import csr_matrix

path = r"D:\manuela\Thesis\results\sn\annotated-muscle-soupxed.h5ad"
f = h5py.File(path, "r")
g = f["raw"]["X"]; shape = tuple(g.attrs["shape"])
X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).tocsc()
genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["raw"]["var"]["_index"][:]])
gidx = {n: int(np.where(genes==n)[0][0]) for n in ["Col3a1","Col4a1","Mmp2","Timp2"]}
print("raw/X max:", round(float(X.max()),2), "(log-norm if <15)")

def cat(col):
    o=f["obs"][col]; cats=np.array([c.decode() if isinstance(c,bytes) else c for c in o["categories"][:]]); return cats[o["codes"][:]]
celltype = cat("celltype")
sample = cat("sample")
f.close()
# sample -> condition (CLAUDE.md: rgzj1=ctrl,2=age,3=sAct,4=DR,5=combi)
smap = {"rgzj1":"ctrl","rgzj2":"age","rgzj3":"sAct","rgzj4":"DR","rgzj5":"combi"}
condition = np.array([smap.get(s, s) for s in sample])

conds = ["ctrl","age","DR","sAct","combi"]
cols = {gn: np.asarray(X[:, gidx[gn]].todense()).ravel() for gn in gidx}
def m_log(ct, gn, c):
    mask=(celltype==ct)&(condition==c); v=cols[gn][mask]; return v.mean(), mask.sum()
def m_logCP(ct, gn, c):
    mask=(celltype==ct)&(condition==c); v=cols[gn][mask]; return np.log1p(np.expm1(v).mean())

ct="FAPs"
print(f"\n== MUSCLE FAPs  (n per cond: " + ", ".join(f"{c}:{int(((celltype==ct)&(condition==c)).sum())}" for c in conds) + ") ==")
for num,den in [("Col3a1","Col4a1"),("Mmp2","Timp2")]:
    print(f"\n  {num}/{den}:")
    for c in conds:
        mn,_=m_log(ct,num,c); md,_=m_log(ct,den,c)
        print(f"    {c:5s}: {num}={mn:.3f} {den}={md:.3f} | ratio(meanLog)={mn/md:.2f}  ratio(log-meanCP)={m_logCP(ct,num,c)/m_logCP(ct,den,c):.2f}")
print("""
THESIS FAPs (4.3): Col3a1/Col4a1 ctrl 2.75, age 1.32, combi 0.50  (Col4a1 0.67 -> 1.25)
                   Mmp2/Timp2 ctrl 1.01, combi 1.76
""")
