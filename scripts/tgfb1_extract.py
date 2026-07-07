# -*- coding: utf-8 -*-
import csv, h5py, numpy as np
from scipy.sparse import csr_matrix

def find_row(path, sym):
    with open(path, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            s = (r.get('Symbol') or r.get('symbol') or r.get('gene_symbol') or '').strip()
            if s == sym: return r
    return None

print("=== BULK Tgfb1 ===")
for tis,age_f,imp_f in [
  ('KIDNEY','results/Kidney_analysis_results/1_DEGs_tables/kidney_res_age_results_with_symbols.csv',
   'results/Kidney_analysis_results/3_Effect_comparison_tables/kidney_intervention_impact_comparison.csv'),
  ('MUSCLE','results/Muscle_v2_analysis_results/1_DEGs_tables/muscle_v2_res_age_results_with_symbols.csv',
   'results/Muscle_v2_analysis_results/3_Effect_comparison_tables/muscle_v2_intervention_impact_comparison.csv')]:
    a=find_row(age_f,'Tgfb1'); im=find_row(imp_f,'Tgfb1')
    la=float(a['log2FoldChange']); pa=float(a['padj'])
    print(f"  {tis}: aging LFC={la:+.3f} padj={pa:.2e} | rescue(combi-vs-age) LFC_Combined={float(im['LFC_Combined']):+.3f}" if im else f"  {tis}: aging only")

print("\n=== snRNA Tgfb1 detection rate (% nuclei >0) ===")
def cats(g):
    c=[x.decode() if isinstance(x,bytes) else x for x in g['categories'][:]]
    return np.array([c[i] if i>=0 else None for i in g['codes'][:]],dtype=object)
CONDS=['ctrl','age','DR','sAct','combi']
def detect(h5,smap,celltypes_match):
    h=h5py.File(h5,'r')
    ct=cats(h['obs']['celltype']); samp=cats(h['obs']['sample'])
    cond=np.array([smap.get(s,s) for s in samp],dtype=object)
    rv=[x.decode() if isinstance(x,bytes) else x for x in h['raw']['var']['_index'][:]]
    rx=h['raw']['X']; M=csr_matrix((rx['data'][:],rx['indices'][:],rx['indptr'][:]),shape=(len(samp),len(rv))).tocsc()
    col=np.asarray(M[:,rv.index('Tgfb1')].todense()).ravel()
    cts=sorted(set(ct))
    target=[c for c in cts if any(k.lower() in str(c).lower() for k in celltypes_match)]
    h.close()
    return col,ct,cond,target,cts
# kidney
col,ct,cond,target,allct=detect('original_data/sn/annotated-aging-soupxed.h5ad',
    {'sg18':'ctrl','sg20':'age','sg24':'sAct','sg26':'DR','sg28':'combi'} , ['IMM','immune'])
print("  kidney celltypes w/ immune:",target, "| (all has IMM?)", [c for c in allct if 'IMM' in str(c) or 'mmune' in str(c)])
for cellt in (target or []):
    vals={cc:round(100*float((col[(ct==cellt)&(cond==cc)]>0).mean()),1) for cc in CONDS}
    print(f"    {cellt}:", [vals[c] for c in CONDS])
# muscle
col,ct,cond,target,allct=detect('original_data/sn/annotated-muscle-soupxed.h5ad',
    {'rgzj1':'ctrl','rgzj2':'age','rgzj3':'sAct','rgzj4':'DR','rgzj5':'combi'}, ['macrophage'])
print("  muscle macrophage label:",target)
for cellt in (target or []):
    vals={cc:round(100*float((col[(ct==cellt)&(cond==cc)]>0).mean()),1) for cc in CONDS}
    print(f"    {cellt}:", [vals[c] for c in CONDS])
