# -*- coding: utf-8 -*-
import pandas as pd, numpy as np, json
from scipy.stats import pearsonr
np.seterr(all='ignore')
KID=r'original_data/bulk/kidney_bulk/deseq2_out'; KIDP=r'original_data/bulk/kidney_bulk/deseq2_out_pairwise'
MUS=r'original_data/bulk/muscle_bulk/v2/deseq2_out'
C=['age','intv','intv1','intv2','combi_ctrl']
def build(dirp):
    a=pd.read_csv(f'{dirp}/res_age_df.csv'); base=a[['gene_id']].copy()
    if 'gene_symbol' in a: base['sym']=a['gene_symbol']
    m=base
    for c in C:
        d=pd.read_csv(f'{dirp}/res_{c}_df.csv')[['gene_id','log2FoldChange','padj']]
        m=m.merge(d.rename(columns={'log2FoldChange':f'lfc_{c}','padj':f'padj_{c}'}),on='gene_id',how='outer')
    return m
def frame(m,thr):
    deg=m[(m['padj_age'].notna())&(m['padj_age']<0.05)&(m['lfc_age'].notna())].copy()
    a=np.sign(deg['lfc_age'])
    r=deg[(np.sign(deg['lfc_intv'])==-a).fillna(False)]
    rtn=r[r['lfc_combi_ctrl'].abs()<thr].copy()
    return deg,r,rtn
def drivers(rtn):
    a=np.sign(rtn['lfc_age'])
    rdr=(np.sign(rtn['lfc_intv1'])==-a).fillna(False); rsa=(np.sign(rtn['lfc_intv2'])==-a).fillna(False)
    n=len(rtn)
    return dict(n=n,dr_excl=int((rdr&~rsa).sum()),sa_excl=int((rsa&~rdr).sum()),
                dual=int((rdr&rsa).sum()),combi_only=int((~rdr&~rsa).sum())), rdr, rsa

def run():
    mk=build(KID); mm=build(MUS)
    kdeg,kr,krtn=frame(mk,0.4)   # kidney 0.4
    mdeg,mr,mrtn=frame(mm,0.5)   # muscle 0.5
    out={}
    kd,krdr,krsa=drivers(krtn); md,mrdr,mrsa=drivers(mrtn)
    out['kidney_drivers@0.4']=kd; out['muscle_drivers@0.5']=md
    # restoration proxy within dual
    def rest(rtn,rdr,rsa):
        dual=rtn[rdr&rsa]; a=np.sign(dual['lfc_age'])
        drs=(-a*dual['lfc_intv1']); sas=(-a*dual['lfc_intv2'])
        return dict(dr_med=round(float(drs.median()),3),sa_med=round(float(sas.median()),3),
                    pct_dr_gt_sa=round(100*float((drs>sas).mean()),1),n_dual=int(len(dual)))
    out['kidney_restore']=rest(krtn,krdr,krsa); out['muscle_restore']=rest(mrtn,mrdr,mrsa)
    # cross-tissue on kidney0.4 set
    kset_sym=set(krtn['sym'].dropna());
    symmap=mk[['gene_id','sym']].dropna().drop_duplicates('gene_id')
    mrtn2=mrtn.merge(symmap,on='gene_id',how='left',suffixes=('','_k'))
    mset_sym=set(mrtn2['sym'].dropna())
    shared_sym=kset_sym & mset_sym
    kset_id=set(krtn['gene_id']); mset_id=set(mrtn['gene_id']); shared_id=kset_id & mset_id
    # consensus + pearson by symbol
    ka=krtn.dropna(subset=['sym']).drop_duplicates('sym').set_index('sym')['lfc_age']
    ma=mrtn2.dropna(subset=['sym']).drop_duplicates('sym').set_index('sym')['lfc_age']
    rows=[(g,ka[g],ma[g]) for g in shared_sym if g in ka.index and g in ma.index]
    sh=pd.DataFrame(rows,columns=['g','k','m']).dropna()
    cons=int((np.sign(sh['k'])==np.sign(sh['m'])).sum()); r,p=pearsonr(sh['k'],sh['m'])
    out['cross@kidney0.4']=dict(shared_by_symbol=len(shared_sym),shared_by_ensembl=len(shared_id),
        pct_of_kidney=round(100*len(shared_sym)/len(kset_sym),1),pct_of_muscle=round(100*len(shared_sym)/len(mset_sym),1),
        consensus=f'{cons}/{len(sh)}',pct_consensus=round(100*cons/len(sh),1),pearson_r=round(float(r),3),pearson_p=float(f'{p:.2e}'))
    return out

runs=[json.dumps(run(),sort_keys=True) for _ in range(3)]
print('3-RUN CONSENSUS:', 'IDENTICAL' if runs[0]==runs[1]==runs[2] else 'MISMATCH')
print(json.dumps(json.loads(runs[0]),indent=1))

# biomarkers ashr vs pairwise (kidney): age LFC + combi_ctrl LFC/padj
print('\n=== §4.9 biomarkers kidney: ashr vs pairwise (MLE) ===')
genes=['Tgfb1','Eda2r','Cdkn1a','Cdkn2a','Gdf15']
for nm,dirp in [('ashr',KID),('pairwise',KIDP)]:
    a=pd.read_csv(f'{dirp}/res_age_df.csv'); cc=pd.read_csv(f'{dirp}/res_combi_ctrl_df.csv')
    for g in genes:
        ra=a[a['gene_symbol']==g] if 'gene_symbol' in a else a.iloc[0:0]
        rc=cc[cc['gene_symbol']==g] if 'gene_symbol' in cc else cc.iloc[0:0]
        if len(ra):
            al=float(ra.iloc[0]['log2FoldChange']); ap=float(ra.iloc[0]['padj'])
            s=f'  {nm:8} {g:8} age_lfc={al:.3f} age_padj={ap:.2e}'
            if len(rc):
                cl=float(rc.iloc[0]['log2FoldChange']); cp=float(rc.iloc[0]['padj'])
                s+=f'  combi_ctrl_lfc={cl:.3f} padj={cp:.2e}'
            print(s)
