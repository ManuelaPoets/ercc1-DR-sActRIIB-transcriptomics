import h5py, numpy as np
from scipy.sparse import csr_matrix

def load(path):
    f = h5py.File(path, "r")
    g = f["raw"]["X"]; shape = tuple(g.attrs["shape"])
    X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).tocsc()
    genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["raw"]["var"]["_index"][:]])
    def cat(c):
        o=f["obs"][c]; cats=np.array([z.decode() if isinstance(z,bytes) else z for z in o["categories"][:]]); return cats[o["codes"][:]]
    celltype=cat("celltype"); sample=cat("sample")
    f.close()
    return X, genes, celltype, sample

conds=["ctrl","age","DR","sAct","combi"]
def det(X,genes,celltype,condition,ct,gn,c):
    if gn not in genes: return None
    gi=int(np.where(genes==gn)[0][0]); v=np.asarray(X[:,gi].todense()).ravel()
    m=(celltype==ct)&(condition==c); return 100*np.mean(v[m]>0) if m.sum() else float("nan")

# KIDNEY
Xk,gk,ctk,sk = load(r"D:\manuela\Thesis\results\sn\annotated-aging-soupxed.h5ad")
condk = sk  # already ctrl/age/...
print("== KIDNEY IMM ==  n:", {c:int(((ctk=='IMM')&(condk==c)).sum()) for c in conds})
for gn in ["H2-Aa","H2-Ab1","H2-Eb1","Adgre1","Csf1r"]:
    print(f"  {gn:7s}: " + "  ".join(f"{c}={det(Xk,gk,ctk,condk,'IMM',gn,c):.1f}%" for c in conds))

# MUSCLE
Xm,gm,ctm,sm = load(r"D:\manuela\Thesis\results\sn\annotated-muscle-soupxed.h5ad")
smap={"rgzj1":"ctrl","rgzj2":"age","rgzj3":"sAct","rgzj4":"DR","rgzj5":"combi"}
condm=np.array([smap.get(s,s) for s in sm])
print("\n== MUSCLE Macrophages ==  n:", {c:int(((ctm=='Macrophages')&(condm==c)).sum()) for c in conds})
for gn in ["Adgre1","Csf1r","H2-Aa","H2-Ab1","H2-Eb1"]:
    print(f"  {gn:7s}: " + "  ".join(f"{c}={det(Xm,gm,ctm,condm,'Macrophages',gn,c):.1f}%" for c in conds))

print("""
THESIS 4.4 KIDNEY IMM: H2-Aa ctrl61.5 age46.1 combi16.0 ; Adgre1 ctrl27.9 sAct37.4 ; Csf1r ctrl32.8 sAct43.4
THESIS 4.4 MUSCLE Macro: Adgre1 57.6/62.2/69.9 ; Csf1r 52.0/60.8/67.5 ; H2-Aa 22.0->30.9 combi9.2 (ctrl/age/combi)
""")
