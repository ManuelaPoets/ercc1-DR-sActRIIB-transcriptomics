import h5py, numpy as np
from scipy.sparse import csr_matrix

path = r"D:\manuela\Thesis\notebooks\aging_interorgan\annotated-interorgan-soupxed.h5ad"
f = h5py.File(path, "r")

# --- raw log-norm matrix (cells x genes) ---
g = f["raw"]["X"]
shape = tuple(g.attrs["shape"])
X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["raw"]["var"]["_index"][:]])
gidx = {name: int(np.where(genes == name)[0][0]) for name in ["Col3a1","Col4a1","Mmp2","Timp2"] if name in genes}
print("genes found:", gidx)
Xc = X.tocsc()

# --- obs categoricals ---
def cat(col):
    o = f["obs"][col]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in o["categories"][:]])
    codes = o["codes"][:]
    return cats[codes]
celltype = cat("celltype")
condition = cat("condition")
f.close()

conds = ["ctrl","age","DR","sAct","combi"]
def mean_expr(ct, gene, cond):
    col = Xc[:, gidx[gene]].toarray().ravel()
    mask = (celltype == ct) & (condition == cond)
    return col[mask].mean() if mask.sum() else float("nan"), int(mask.sum())

def report(ct, label):
    print(f"\n==== {label}  (celltype == '{ct}') ====")
    n = {c: int(((celltype==ct)&(condition==c)).sum()) for c in conds}
    print("n cells per condition:", n)
    for num, den in [("Col3a1","Col4a1"), ("Mmp2","Timp2")]:
        print(f"\n  {num}/{den} ratio (mean log-norm expression):")
        for c in conds:
            mn,_ = mean_expr(ct, num, c); md,_ = mean_expr(ct, den, c)
            ratio = mn/md if md else float("nan")
            print(f"    {c:5s}: {num}={mn:.3f}  {den}={md:.3f}  ratio={ratio:.2f}")

report("FIB", "KIDNEY fibroblasts")
report("FAPs", "MUSCLE FAPs")

print("""
THESIS (4.3) — KIDNEY FIB: Col3a1/Col4a1 ctrl 2.90, age 3.31, DR 3.15, sAct 2.89, combi 1.16
   (Col4a1 0.80->0.73; Col3a1 2.31->2.41 age; combi Col3a1 0.87, Col4a1 0.75; DR Col3a1 1.60)
   Mmp2/Timp2 ctrl 0.23, age 1.41 (Mmp2 0.272->1.213; Timp2 1.17->0.86), combi 0.83 (Mmp2 0.52, Timp2 0.62)
THESIS — MUSCLE FAPs: Col3a1/Col4a1 ctrl 2.75, age 1.32 (Col4a1 0.67->1.25), combi 0.50
   Mmp2/Timp2 combi 1.76 vs ctrl 1.01
""")
