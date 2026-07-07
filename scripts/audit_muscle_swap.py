# -*- coding: utf-8 -*-
"""Full audit of muscle DR/sAct snRNA values. CORRECT map: rgzj3=sAct, rgzj4=DR."""
import h5py, numpy as np, json
from scipy.sparse import csr_matrix
MAP={'rgzj1':'ctrl','rgzj2':'age','rgzj3':'sAct','rgzj4':'DR','rgzj5':'combi'}  # CORRECT
CONDS=['ctrl','age','DR','sAct','combi']
def cats(g):
    c=[x.decode() if isinstance(x,bytes) else x for x in g['categories'][:]]
    return np.array([c[i] if i>=0 else None for i in g['codes'][:]],dtype=object)
h=h5py.File(r'original_data/sn/annotated-muscle-soupxed.h5ad','r')
ct=cats(h['obs']['celltype']); rg=cats(h['obs']['sample'])
cond=np.array([MAP[s] for s in rg],dtype=object)
rv=[x.decode() if isinstance(x,bytes) else x for x in h['raw']['var']['_index'][:]]
rx=h['raw']['X']; M=csr_matrix((rx['data'][:],rx['indices'][:],rx['indptr'][:]),shape=(len(rg),len(rv))).tocsc()
def colv(g): return np.asarray(M[:,rv.index(g)].todense()).ravel() if g in rv else None
def mean_by(g,cellt):
    c=colv(g); return {cc:round(float(c[(ct==cellt)&(cond==cc)].mean()),3) for cc in CONDS}
def det_by(g,cellt):
    c=colv(g); return {cc:round(100*float((c[(ct==cellt)&(cond==cc)]>0).mean()),2) for cc in CONDS}

print('=== CORRECT muscle values (rgzj3=sAct, rgzj4=DR) ===')
# D6 FAP Col3a1/Col4a1 ratio
c3=mean_by('Col3a1','FAPs'); c4=mean_by('Col4a1','FAPs')
ratio={cc:round(c3[cc]/c4[cc],2) for cc in CONDS}
print('D6 FAP Col3a1/Col4a1 ratio [ctrl,age,DR,sAct,combi]:', [ratio[c] for c in CONDS], '  (figure has [2.14,1.32,1.37,0.86,0.68])')
# D6 FAP Mmp2/Timp2
m2=mean_by('Mmp2','FAPs'); t2=mean_by('Timp2','FAPs')
mr={cc:round(m2[cc]/t2[cc],2) for cc in CONDS}
print('D6 FAP Mmp2/Timp2 ratio   [ctrl,age,DR,sAct,combi]:', [mr[c] for c in CONDS], '  (figure has [0.84,1.31,1.23,1.08,1.23])')
# D8 MHC-II muscle macrophages
print('D8 Macrophage detection [ctrl,age,DR,sAct,combi]:')
for g in ['H2-Aa','H2-Ab1','H2-Eb1','Adgre1','Csf1r']:
    d=det_by(g,'Macrophages'); print(f'   {g:7}:', [d[c] for c in CONDS])
h.close()
