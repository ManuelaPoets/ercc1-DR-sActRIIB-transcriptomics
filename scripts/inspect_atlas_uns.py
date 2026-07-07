import h5py

def show(f, path):
    try:
        g = f[path]
    except KeyError:
        print(path, "MISSING"); return
    print("==", path, "==")
    for k, v in g.attrs.items():
        print("  attr", k, "=", v)
    if isinstance(g, h5py.Group):
        for k in g.keys():
            node = g[k]
            if isinstance(node, h5py.Dataset):
                print("   ", k, "=", node[()])
            else:
                print("   ", k, "(group):")
                for kk in node.keys():
                    nn = node[kk]
                    if isinstance(nn, h5py.Dataset):
                        print("       ", kk, "=", nn[()])
                    for ak, av in nn.attrs.items():
                        print("        attr", kk, ak, "=", av)

import sys
path = sys.argv[1] if len(sys.argv) > 1 else r"D:\manuela\Thesis\results\sn\annotated-aging-soupxed.h5ad"
f = h5py.File(path, "r")
print("FILE:", path)
print("UNS keys:", list(f["uns"].keys()))
for p in ["uns/neighbors", "uns/leiden", "uns/pca", "uns/hvg", "uns/log1p", "uns/umap"]:
    show(f, p)
print("\nOBSM keys:", list(f["obsm"].keys()))
print("X_pca_harmony present:", "X_pca_harmony" in f["obsm"])
# n cells
print("n_obs (X rows attr):", dict(f["X"].attrs) if "X" in f else "no X")
f.close()
