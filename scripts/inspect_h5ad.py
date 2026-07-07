import h5py, numpy as np, sys

path = sys.argv[1] if len(sys.argv) > 1 else r"D:\manuela\Thesis\notebooks\aging_interorgan\annotated-interorgan-soupxed.h5ad"
f = h5py.File(path, "r")
print("FILE:", path)
print("top-level keys:", list(f.keys()))

def shape_of(g):
    if isinstance(g, h5py.Dataset): return ("dense", g.shape, g.dtype)
    if isinstance(g, h5py.Group):
        if "data" in g and "indptr" in g:
            enc = g.attrs.get("encoding-type", b"?")
            return ("sparse", g.attrs.get("shape", "?"), enc)
        return ("group", list(g.keys()))
    return "?"
print("X:", shape_of(f["X"]))
if "raw" in f:
    print("raw keys:", list(f["raw"].keys()))
    print("raw/X:", shape_of(f["raw"]["X"]))

def obs_cols(f):
    o = f["obs"]
    print("\nobs keys:", list(o.keys()))
    for k in o.keys():
        item = o[k]
        if isinstance(item, h5py.Group) and "categories" in item:  # categorical
            cats = item["categories"][:]
            cats = [c.decode() if isinstance(c, bytes) else c for c in cats]
            print(f"  [cat] {k}: {len(cats)} cats -> {cats[:30]}")
        elif isinstance(item, h5py.Dataset):
            print(f"  [arr] {k}: shape={item.shape} dtype={item.dtype}")
obs_cols(f)

# var / raw-var names location
for grp in ["var", "raw"]:
    if grp == "raw" and "raw" in f and "var" in f["raw"]:
        v = f["raw"]["var"]; label="raw/var"
    elif grp == "var":
        v = f["var"]; label="var"
    else:
        continue
    idxname = v.attrs.get("_index", b"_index")
    idxname = idxname.decode() if isinstance(idxname, bytes) else idxname
    print(f"\n{label} _index='{idxname}', keys={list(v.keys())}")
f.close()
