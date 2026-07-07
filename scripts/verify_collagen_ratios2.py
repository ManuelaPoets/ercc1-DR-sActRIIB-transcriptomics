import h5py, numpy as np
from scipy.sparse import csr_matrix

path = r"D:\manuela\Thesis\notebooks\aging_interorgan\annotated-interorgan-soupxed.h5ad"
f = h5py.File(path, "r")
g = f["raw"]["X"]; shape = tuple(g.attrs["shape"])
X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).tocsc()
genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["raw"]["var"]["_index"][:]])
gidx = {n: int(np.where(genes==n)[0][0]) for n in ["Col3a1","Col4a1","Mmp2","Timp2"]}
def cat(col):
    o=f["obs"][col]; cats=np.array([c.decode() if isinstance(c,bytes) else c for c in o["categories"][:]]); return cats[o["codes"][:]]
celltype=cat("celltype"); condition=cat("condition")
# also check max value of raw to see if it's log (<~10) or counts (large)
print("raw/X max value:", X.max(), "| looks like", "log-norm" if X.max()<15 else "counts")
f.close()

conds=["ctrl","age","DR","sAct","combi"]
def col(gene): return np.asarray(X[:, gidx[gene]].todense()).ravel()
cols={g_:col(g_) for g_ in gidx}

def stats(ct, gene, cond):
    m=(celltype==ct)&(condition==cond); v=cols[gene][m]
    return dict(mean_log=v.mean(), log_of_meanCP=np.log1p(np.expm1(v).mean()), meanCP=np.expm1(v).mean())

def report(ct,label):
    print(f"\n==== {label} ({ct}) ====")
    for num,den in [("Col3a1","Col4a1"),("Mmp2","Timp2")]:
        print(f"  {num}/{den}:")
        print("    cond   |  meanLog  log(meanCP) |  ratio_meanLog  ratio_log(meanCP)")
        for c in conds:
            sn=stats(ct,num,c); sd=stats(ct,den,c)
            print(f"    {c:5s} | {num}:{sn['log_of_meanCP']:.3f} {den}:{sd['log_of_meanCP']:.3f} | "
                  f"{sn['mean_log']/sd['mean_log']:.2f}  {sn['log_of_meanCP']/sd['log_of_meanCP']:.2f}")

report("FIB","KIDNEY fibroblasts")
report("FAPs","MUSCLE FAPs")
print("\nTHESIS FIB Col3a1/Col4a1: ctrl2.90 age3.31 DR3.15 sAct2.89 combi1.16 (Col4a1 0.80/0.73; Col3a1 2.31/2.41)")
print("THESIS FIB Mmp2/Timp2: ctrl0.23 age1.41 combi0.83 (Mmp2 0.272/1.213; Timp2 1.17/0.86)")
print("THESIS FAPs Col3a1/Col4a1: ctrl2.75 age1.32 combi0.50 ; Mmp2/Timp2 ctrl1.01 combi1.76")
