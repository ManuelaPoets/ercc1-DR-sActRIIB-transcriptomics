# -*- coding: utf-8 -*-
import json, re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
nb = json.load(open(r"D:\manuela\Thesis\notebooks\results_plots_clean.ipynb", encoding="utf-8"))
cells = nb["cells"]

SN_ATLAS = re.compile(r'muscle_umap_coordinates|muscle_celltype_condition')
OLD = re.compile(r'Muscle_old')
V2  = re.compile(r'Muscle_v2|MUSCLE_V2')
FIGID = re.compile(r"(?:FILENAME\w*\s*=\s*['\"]|save_fig\([^,]+,\s*['\"]|f?['\"])(Fig[_\. ]?3[_\.]?\w+|Table[_\. ]?3\w+|S3\.\w+|S4\.\w+)", re.I)
HDR = re.compile(r'^#+\s*(Fig|Table|Supp|S3|S4|#?\s*3\.\d)', re.I)

last_md = ""
rows = []
for i, c in enumerate(cells):
    src = "".join(c.get("source", []))
    if c["cell_type"] == "markdown":
        # capture figure-ish headers
        first = src.strip().splitlines()[0] if src.strip() else ""
        last_md = first
        continue
    if c["cell_type"] != "code" or not src.strip():
        continue
    has_old = bool(OLD.search(src))
    has_v2  = bool(V2.search(src))
    has_sn_atlas = bool(SN_ATLAS.search(src))
    # bulk/hybrid old = old reference that is NOT only the sn-atlas files
    old_lines = [l for l in src.splitlines() if 'Muscle_old' in l]
    old_bulk = any(not SN_ATLAS.search(l) for l in old_lines)
    old_sn_only = len(old_lines) > 0 and all(SN_ATLAS.search(l) for l in old_lines)
    # classify
    if has_old and has_v2:
        cls = "MIXED(old+v2)"
    elif has_v2:
        cls = "V2"
    elif old_sn_only:
        cls = "OLD-SN-ATLAS"
    elif old_bulk:
        cls = "OLD-BULK/HYBRID"
    elif has_old:
        cls = "OLD-other"
    else:
        cls = "no-muscle"
    figids = sorted(set(m.group(1) for m in FIGID.finditer(src)))
    rows.append((i, cls, last_md[:48], ";".join(figids)[:60]))

print(f"{'cell':>4} {'class':<16} {'preceding header':<50} figIDs")
print("-"*120)
for r in rows:
    print(f"{r[0]:>4} {r[1]:<16} {r[2]:<50} {r[3]}")

print("\n===== SUMMARY =====")
from collections import Counter
cnt = Counter(r[1] for r in rows)
for k,v in cnt.items(): print(f"  {k}: {v}")
print("\nOLD-BULK/HYBRID cells (drop candidates):", [r[0] for r in rows if r[1]=="OLD-BULK/HYBRID"])
print("OLD-SN-ATLAS cells (keep as-is):", [r[0] for r in rows if r[1]=="OLD-SN-ATLAS"])
print("MIXED cells (inspect):", [r[0] for r in rows if r[1].startswith("MIXED")])
