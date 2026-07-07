# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
np.seterr(all='ignore')
KID=r'results/kidney_bulk/deseq2_out'; KIDP=r'results/kidney_bulk/deseq2_out_pairwise'
MUS=r'results/muscle_bulk/v2/deseq2_out'
C=['age','intv','intv1','intv2','combi_ctrl']
def build(dirp):
    a=pd.read_csv(f'{dirp}/res_age_df.csv')
    base=a[['gene_id']].copy()
    if 'gene_symbol' in a: base['sym']=a['gene_symbol']
    m=base
    for c in C:
        d=pd.read_csv(f'{dirp}/res_{c}_df.csv')[['gene_id','log2FoldChange','padj']]
        m=m.merge(d.rename(columns={'log2FoldChange':f'lfc_{c}','padj':f'padj_{c}'}),on='gene_id',how='outer')
    return m
def frame(m,thr):
    deg=m[(m['padj_age'].notna())&(m['padj_age']<0.05)&(m['lfc_age'].notna())].copy()
    a=np.sign(deg['lfc_age'])
    rescued=(np.sign(deg['lfc_intv'])==-a).fillna(False)
    r=deg[rescued]
    rtn=r[r['lfc_combi_ctrl'].abs()<thr].copy()
    return deg,r,rtn
def drivers(rtn):
    a=np.sign(rtn['lfc_age'])
    rdr=(np.sign(rtn['lfc_intv1'])==-a).fillna(False)
    rsa=(np.sign(rtn['lfc_intv2'])==-a).fillna(False)
    n=len(rtn)
    d=dict(dr_excl=int((rdr&~rsa).sum()),sa_excl=int((rsa&~rdr).sum()),
           dual=int((rdr&rsa).sum()),combi_only=int((~rdr&~rsa).sum()))
    d.update({k+'_pct':round(100*v/n,1) for k,v in list(d.items())})
    d['n']=n
    return d,rdr,rsa

print('### KIDNEY drivers under thr=0.4 (should match thesis 546/108/87/1220):')
mk=build(KID); kdeg,kr,krtn=frame(mk,0.4)
kd,krdr,krsa=drivers(krtn); print('  ',kd)
print('### MUSCLE drivers under thr=0.5 (thesis 211/20/14/246):')
mm=build(MUS); mdeg,mr,mrtn=frame(mm,0.5)
md,mrdr,mrsa=drivers(mrtn); print('  ',md)

# Restoration-score proxy within DUAL-rescued genes: reversal magnitude in aging-opposite direction
def restore_scores(rtn,rdr,rsa):
    dual=rtn[rdr&rsa].copy()
    a=np.sign(dual['lfc_age'])
    dr_s=(-a*dual['lfc_intv1'])  # how far DR moves opposite to aging
    sa_s=(-a*dual['lfc_intv2'])
    return dual,dr_s,sa_s
for name,rtn,rdr,rsa in [('KIDNEY',krtn,krdr,krsa),('MUSCLE',mrtn,mrdr,mrsa)]:
    dual,dr_s,sa_s=restore_scores(rtn,rdr,rsa)
    pct=round(100*(dr_s>sa_s).mean(),1)
    print(f'### {name} restoration proxy (-sign(age)*lfc): DR median={dr_s.median():.3f} sAct median={sa_s.median():.3f} DR>sAct in {pct}% of {len(dual)} dual')

# Kidney Tgfb1 in ashr vs pairwise
for nm,dirp in [('ashr',KID),('pairwise',KIDP)]:
    a=pd.read_csv(f'{dirp}/res_age_df.csv')
    row=a[a['gene_symbol']=='Tgfb1'] if 'gene_symbol' in a else None
    if row is not None and len(row):
        print(f'Kidney Tgfb1 {nm}: age lfc={row.iloc[0][\"log2FoldChange\"]:.3f} padj={row.iloc[0][\"padj\"]:.2e}')
