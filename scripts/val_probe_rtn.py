# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
np.seterr(all='ignore')
KID=r'results/kidney_bulk/deseq2_out'; MUS=r'results/muscle_bulk/v2/deseq2_out'
C=['age','intv','intv1','intv2','combi_ctrl']
def build(dirp):
    base=pd.read_csv(f'{dirp}/res_age_df.csv')[['gene_id']].copy()
    m=base
    for c in C:
        d=pd.read_csv(f'{dirp}/res_{c}_df.csv')[['gene_id','log2FoldChange','padj']]
        m=m.merge(d.rename(columns={'log2FoldChange':f'lfc_{c}','padj':f'padj_{c}'}),on='gene_id',how='outer')
    return m
def stage1(m):
    deg=m[(m['padj_age'].notna())&(m['padj_age']<0.05)&(m['lfc_age'].notna())].copy()
    a=np.sign(deg['lfc_age'])
    rescued=(np.sign(deg['lfc_intv'])==-a).fillna(False)
    return deg,deg[rescued]
for name,dirp,target in [('KIDNEY',KID,1961),('MUSCLE',MUS,491)]:
    m=build(dirp); deg,r=stage1(m)
    print(f'=== {name}  aging DEG={len(deg)}  stage1 rescued={len(r)}  target RTN={target}')
    cc=r['lfc_combi_ctrl']; ccp=r['padj_combi_ctrl']
    for thr in [0.3,0.4,0.5,0.585,0.6,1.0]:
        print(f'   |combi_ctrl_lfc|<{thr}: {int((cc.abs()<thr).sum())}')
    print(f'   combi_ctrl padj>=0.05 (ns vs ctrl): {int((ccp.isna()|(ccp>=0.05)).sum())}')
    print(f'   |combi_ctrl|<0.5 AND padj>=0.05: {int(((cc.abs()<0.5)&(ccp.isna()|(ccp>=0.05))).sum())}')
    print(f'   |combi_ctrl|<1.0 AND padj>=0.05: {int(((cc.abs()<1.0)&(ccp.isna()|(ccp>=0.05))).sum())}')
    # |combi_ctrl| < |age|/2 ?  relative restoration
    rel=(cc.abs() < 0.5*r['lfc_age'].abs())
    print(f'   |combi_ctrl| < 0.5*|age_lfc|: {int(rel.sum())}')
