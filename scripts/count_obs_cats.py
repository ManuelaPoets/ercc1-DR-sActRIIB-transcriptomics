import h5py, sys

def cats(f, col):
    try:
        node = f[f"obs/{col}"]
    except KeyError:
        return f"{col}: MISSING"
    if isinstance(node, h5py.Group):  # categorical
        c = node["categories"][()]
        c = [x.decode() if isinstance(x, bytes) else x for x in c]
        return f"{col}: {len(c)} cats -> {c}"
    else:
        import numpy as np
        v = node[()]
        u = sorted(set(v.tolist()))
        return f"{col}: {len(u)} unique -> {u[:40]}"

path = sys.argv[1]
f = h5py.File(path, "r")
print("FILE:", path)
print("obs keys:", list(f["obs"].keys()))
for col in ["leiden", "clusters", "subclusters", "mergedClusters", "celltype", "subtype", "supertype"]:
    if col in f["obs"]:
        print(" ", cats(f, col))
f.close()
