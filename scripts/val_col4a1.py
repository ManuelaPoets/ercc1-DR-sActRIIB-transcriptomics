# -*- coding: utf-8 -*-
"""Confirm §4.3 Col4a1 (and coupled Col3a1 / ratio) values.
Metric: mean of raw/X log-norm within celltype x condition. Run 3x for consensus."""
import h5py, numpy as np
from scipy.sparse import csr_matrix

def decode(arr):
    return [x.decode() if isinstance(x,bytes) else x for x in arr]

def load_group_means(path, celltype, genes, sample_map=None):
    h=h5py.File(path,'r')
    # obs categorical
    def cat_col(name):
        g=h['obs'][name]
        codes=g['codes'][:]; cats=decode(g['categories'][:])
        return np.array([cats[c] if c>=0 else None for c in codes],dtype=object)
    ct=cat_col('celltype'); samp=cat_col('sample')
    if sample_map: samp=np.array([sample_map.get(s,s) for s in samp],dtype=object)
    # raw var index
    rv=decode(h['raw']['var']['_index'][:])
    gidx={g:rv.index(g) for g in genes if g in rv}
    # build CSR raw/X
    rx=h['raw']['X']
    data=rx['data'][:]; indices=rx['indices'][:]; indptr=rx['indptr'][:]
    ncells=len(indptr)-1; ngenes=len(rv)
    M=csr_matrix((data,indices,indptr),shape=(ncells,ngenes))
    h.close()
    conds=['ctrl','age','DR','sAct','combi']
    out={}
    for g in genes:
        if g not in gidx: out[g]={c:None for c in conds}; out[g]['_missing']=True; continue
        col=np.asarray(M[:,gidx[g]].todense()).ravel()
        res={}
        for c in conds:
            mask=(ct==celltype)&(samp==c)
            res[c]=round(float(col[mask].mean()),3) if mask.sum()>0 else None
            res[c+'_n']=int(mask.sum())
        out[g]=res
    return out

MUS_MAP={'rgzj1':'ctrl','rgzj2':'age','rgzj3':'sAct','rgzj4':'DR','rgzj5':'combi'}
KPATH=r'original_data/sn/annotated-aging-soupxed.h5ad'
MPATH=r'original_data/sn/annotated-muscle-soupxed.h5ad'
def report():
    print('=== KIDNEY FIB ===')
    k=load_group_means(KPATH,'FIB',['Col4a1','Col3a1'])
    print('  cells/cond:',{c:k['Col4a1'][c+'_n'] for c in ['ctrl','age','DR','sAct','combi']})
    for g in ['Col4a1','Col3a1']:
        print(f'  {g}:',{c:k[g][c] for c in ['ctrl','age','DR','sAct','combi']})
    r={c:(round(k['Col3a1'][c]/k['Col4a1'][c],3) if k['Col4a1'][c] else None) for c in ['ctrl','age','DR','sAct','combi']}
    print('  Col3a1/Col4a1 ratio:',r)
    print('  Col4a1 ctrl->age %:', round(100*(k['Col4a1']['age']-k['Col4a1']['ctrl'])/k['Col4a1']['ctrl'],1))
    print('=== MUSCLE FAPs ===')
    m=load_group_means(MPATH,'FAPs',['Col4a1','Col3a1'],MUS_MAP)
    print('  cells/cond:',{c:m['Col4a1'][c+'_n'] for c in ['ctrl','age','DR','sAct','combi']})
    for g in ['Col4a1','Col3a1']:
        print(f'  {g}:',{c:m[g][c] for c in ['ctrl','age','DR','sAct','combi']})
    rm={c:(round(m['Col3a1'][c]/m['Col4a1'][c],3) if m['Col4a1'][c] else None) for c in ['ctrl','age','DR','sAct','combi']}
    print('  Col3a1/Col4a1 ratio:',rm)
    return k,m

import json
runs=[]
for i in range(3):
    k,m=report() if i==0 else (None,None)
    # for consensus, recompute silently runs 2,3
    if i>0:
        k=load_group_means(KPATH,'FIB',['Col4a1','Col3a1'])
        m=load_group_means(MPATH,'FAPs',['Col4a1','Col3a1'],MUS_MAP)
    runs.append(json.dumps([k,m],sort_keys=True,default=str))
print('=== 3-RUN CONSENSUS:', 'ALL IDENTICAL' if runs[0]==runs[1]==runs[2] else 'MISMATCH')
