# -*- coding: utf-8 -*-
"""Validate all snRNA Discussion claims (§4.3-4.8): detection rates (% nuclei>0),
mean log-norm (raw/X), and composition. Run 3x for consensus."""
import h5py, numpy as np, json
from scipy.sparse import csr_matrix
KP=r'original_data/sn/annotated-aging-soupxed.h5ad'
MP=r'original_data/sn/annotated-muscle-soupxed.h5ad'
MUS_MAP={'rgzj1':'ctrl','rgzj2':'age','rgzj3':'sAct','rgzj4':'DR','rgzj5':'combi'}
CONDS=['ctrl','age','DR','sAct','combi']
def dec(a): return [x.decode() if isinstance(x,bytes) else x for x in a]
def load(path, sample_map=None):
    h=h5py.File(path,'r')
    def catcol(n):
        g=h['obs'][n]; codes=g['codes'][:]; cats=dec(g['categories'][:])
        return np.array([cats[c] if c>=0 else None for c in codes],dtype=object)
    ct=catcol('celltype'); samp=catcol('sample')
    if sample_map: samp=np.array([sample_map.get(s,s) for s in samp],dtype=object)
    rv=dec(h['raw']['var']['_index'][:])
    rx=h['raw']['X']; M=csr_matrix((rx['data'][:],rx['indices'][:],rx['indptr'][:]),shape=(len(samp),len(rv)))
    h.close()
    return ct,samp,rv,M.tocsc()
def col(M,rv,g):
    if g not in rv: return None
    return np.asarray(M[:,rv.index(g)].todense()).ravel()
def detect(ct,samp,M,rv,celltype,gene):  # % nuclei >0
    c=col(M,rv,gene)
    if c is None: return {x:'NA' for x in CONDS}
    out={}
    for cond in CONDS:
        m=(ct==celltype)&(samp==cond); out[cond]=round(100*float((c[m]>0).mean()),1) if m.sum() else None
    return out
def meanln(ct,samp,M,rv,celltype,gene):  # mean log-norm
    c=col(M,rv,gene)
    if c is None: return {x:'NA' for x in CONDS}
    return {cond:(round(float(c[(ct==celltype)&(samp==cond)].mean()),3) if ((ct==celltype)&(samp==cond)).sum() else None) for cond in CONDS}
def ncells(ct,samp,celltype):
    return {cond:int(((ct==celltype)&(samp==cond)).sum()) for cond in CONDS}
def comp_pct(ct,samp,celltype):  # % of nuclei in this celltype per condition
    out={}
    for cond in CONDS:
        tot=(samp==cond).sum(); out[cond]=round(100*((ct==celltype)&(samp==cond)).sum()/tot,1) if tot else None
    return out

def run():
    ckt,csamp,crv,CM=load(KP)
    mkt,msamp,mrv,MM=load(MP, MUS_MAP)
    R={}
    # --- §4.3 kidney FIB Mmp2/Timp2 (mean), Lox (detect), Col1a1 (mean) ---
    R['K_FIB_Mmp2_mean']=meanln(ckt,csamp,CM,crv,'FIB','Mmp2')
    R['K_FIB_Timp2_mean']=meanln(ckt,csamp,CM,crv,'FIB','Timp2')
    mm=R['K_FIB_Mmp2_mean']; tt=R['K_FIB_Timp2_mean']
    R['K_FIB_Mmp2_Timp2_ratio']={c:(round(mm[c]/tt[c],3) if tt[c] else None) for c in CONDS}
    R['K_FIB_Lox_detect']=detect(ckt,csamp,CM,crv,'FIB','Lox')
    R['K_FIB_Col1a1_mean']=meanln(ckt,csamp,CM,crv,'FIB','Col1a1')
    # --- §4.3/4.6 muscle FAP Mmp2/Timp2 ratio ---
    fm=meanln(mkt,msamp,MM,mrv,'FAPs','Mmp2'); ft=meanln(mkt,msamp,MM,mrv,'FAPs','Timp2')
    R['M_FAP_Mmp2_mean']=fm; R['M_FAP_Timp2_mean']=ft
    R['M_FAP_Mmp2_Timp2_ratio']={c:(round(fm[c]/ft[c],3) if ft[c] else None) for c in CONDS}
    # --- §4.4 kidney IMM MHC-II + macrophage markers (detection) ---
    for g in ['H2-Aa','H2-Ab1','H2-Eb1','Adgre1','Csf1r']:
        R[f'K_IMM_{g}_detect']=detect(ckt,csamp,CM,crv,'IMM',g)
    # --- §4.4 muscle Macrophages markers (detection) ---
    for g in ['Adgre1','Csf1r','H2-Aa']:
        R[f'M_Mac_{g}_detect']=detect(mkt,msamp,MM,mrv,'Macrophages',g)
    # --- §4.5 PT-1 vs PT-2 segment markers (detection; thesis quotes PT-1 vs PT-2, pooled across conditions) ---
    seg=['Slc7a13','Slc22a13','Cyp7b1','Atp11a','Aqp1','Slc34a1','Slc5a2','Slc5a12','Havcr1','Cdkn1a','Vcam1','Krt20']
    for g in seg:
        c=col(CM,crv,g)
        if c is None: R[f'K_PT_{g}_overall']={'PT-1':'NA','PT-2':'NA'}; continue
        R[f'K_PT_{g}_overall']={ctp:round(100*float((c[ckt==ctp]>0).mean()),1) for ctp in ['PT-1','PT-2']}
    # --- §4.5 PT-2 & PT-1 Vcam1/Havcr1/Egfr by condition (detection) ---
    R['K_PT2_Vcam1_detect']=detect(ckt,csamp,CM,crv,'PT-2','Vcam1')
    R['K_PT1_Vcam1_detect']=detect(ckt,csamp,CM,crv,'PT-1','Vcam1')
    R['K_PT2_Havcr1_detect']=detect(ckt,csamp,CM,crv,'PT-2','Havcr1')
    R['K_PT1_Egfr_detect']=detect(ckt,csamp,CM,crv,'PT-1','Egfr')
    R['K_PT1_Egfr_mean']=meanln(ckt,csamp,CM,crv,'PT-1','Egfr')
    # --- §4.4 Tek in glomerular EC by condition (detection) ---
    R['K_gEC_Tek_detect']=detect(ckt,csamp,CM,crv,'EC-1(gEC)','Tek')
    # --- §4.4/4.3 Pdgfra/Pdgfrb localisation (overall detection across celltypes, find top) ---
    for g in ['Pdgfra','Pdgfrb']:
        c=col(CM,crv,g)
        if c is None: continue
        cts=sorted(set(ckt[ckt!=None]))
        dd={ctp:round(100*float((c[ckt==ctp]>0).mean()),1) for ctp in cts}
        R[f'K_{g}_top4_detect']=dict(sorted(dd.items(),key=lambda kv:-kv[1])[:4])
    # --- §4.5/§4.8 composition ---
    R['K_PT2_ncells']=ncells(ckt,csamp,'PT-2')
    R['K_PT1_ncells']=ncells(ckt,csamp,'PT-1')
    R['K_PT1_pct']=comp_pct(ckt,csamp,'PT-1')
    R['K_FIB_pct']=comp_pct(ckt,csamp,'FIB')
    R['K_vSMCMC_pct']=comp_pct(ckt,csamp,'vSMC/MC')
    return R

runs=[json.dumps(run(),sort_keys=True,default=str) for _ in range(3)]
print('3-RUN CONSENSUS:', 'IDENTICAL' if runs[0]==runs[1]==runs[2] else 'MISMATCH')
print(json.dumps(json.loads(runs[0]),indent=1))
