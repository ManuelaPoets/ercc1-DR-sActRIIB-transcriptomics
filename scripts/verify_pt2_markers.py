import h5py, numpy as np
from scipy.sparse import csr_matrix

path = r"D:\manuela\Thesis\results\sn\annotated-aging-soupxed.h5ad"
f = h5py.File(path, "r")
g = f["raw"]["X"]; shape = tuple(g.attrs["shape"])
X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).tocsc()
genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["raw"]["var"]["_index"][:]])
need = ["Havcr1","Cdkn1a","Slc34a1","Aqp1","Vcam1","Egfr"]
gidx = {n:int(np.where(genes==n)[0][0]) for n in need if n in genes}
print("genes found:", list(gidx))
def cat(c):
    o=f["obs"][c]; cats=np.array([z.decode() if isinstance(z,bytes) else z for z in o["categories"][:]]); return cats[o["codes"][:]]
celltype=cat("celltype"); condition=cat("sample")
f.close()
cols={n:np.asarray(X[:,gidx[n]].todense()).ravel() for n in gidx}
conds=["ctrl","age","DR","sAct","combi"]

def det(ct,gn,c):
    m=(celltype==ct)&(condition==c); v=cols[gn][m]; return 100*np.mean(v>0) if m.sum() else float("nan")
def meanlog(ct,gn,c):
    m=(celltype==ct)&(condition==c); v=cols[gn][m]; return v.mean() if m.sum() else float("nan")

print("\n== PT-2 cell counts per condition (thesis: ctrl738 age186 sAct49 combi185) ==")
print({c:int(((celltype=='PT-2')&(condition==c)).sum()) for c in conds})
print("== PT-1 cell counts ==")
print({c:int(((celltype=='PT-1')&(condition==c)).sum()) for c in conds})

print("\n== detection rate (% nonzero) ==")
for ct,gn in [("PT-2","Havcr1"),("PT-2","Cdkn1a"),("PT-2","Slc34a1"),("PT-1","Slc34a1"),
              ("PT-2","Aqp1"),("PT-1","Aqp1"),("PT-1","Vcam1"),("PT-2","Vcam1")]:
    print(f"  {gn:8s} in {ct}: " + "  ".join(f"{c}={det(ct,gn,c):.1f}%" for c in conds))

print("\n== Egfr mean log-norm in PT-1 (thesis ctrl0.46 -> age0.94) ==")
print("  " + "  ".join(f"{c}={meanlog('PT-1','Egfr',c):.2f}" for c in conds))

print("""
THESIS 4.5: Havcr1 PT-2 3.8/9.1/-/-/20.0 ; Cdkn1a PT-2 0.4/17.2/-/-/14.1 ; Slc34a1 PT-2 51.9/44.1/-/-/13.5 (PT-1 ctrl 100)
            Aqp1 PT-2 73.8 vs PT-1 52.4 ; Vcam1 PT-1 0.8->8.1 combi3.8 ; PT-2 1.4->9.1 combi0.5
""")
