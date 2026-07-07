import h5py, numpy as np
from scipy.sparse import csr_matrix

# Seurat v4 FindMarkers avg_log2FC: log2( (mean(expm1(lognorm[g1]))+1) / (mean(expm1(lognorm[g2]))+1) )
path = r"D:\manuela\Thesis\results\sn\annotated-aging-soupxed.h5ad"
f = h5py.File(path, "r")
g = f["raw"]["X"]; shape = tuple(g.attrs["shape"])
X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).tocsc()
genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["raw"]["var"]["_index"][:]])
o=f["obs"]["celltype"]; cats=np.array([z.decode() if isinstance(z,bytes) else z for z in o["categories"][:]]); celltype=cats[o["codes"][:]]
f.close()

def avglog2fc(gene, ct):
    gi = np.where(genes==gene)[0]
    if not len(gi): return None
    v = np.asarray(X[:, int(gi[0])].todense()).ravel()
    e = np.expm1(v)
    m1 = e[celltype==ct].mean(); m2 = e[celltype!=ct].mean()
    return np.log2((m1+1)/(m2+1))

print("Seurat-style avg_log2FC (cell type vs rest), kidney atlas:")
for gene, ct, thesis in [("Flt1","EC-2",2.57),("Flt1","EC-3",4.01),("Gstm1","PT-2",3.98),("Gstm1","FIB",5.16)]:
    val = avglog2fc(gene, ct)
    print(f"  {gene} in {ct}: computed={val:.2f}  thesis={thesis}")
