# -*- coding: utf-8 -*-
"""Validate Discussion bulk RNA-seq claims by independent reproduction from res_*_df.csv.
Each check is run 3x; we assert all 3 runs agree (consensus)."""
import pandas as pd, numpy as np, json
np.seterr(all='ignore')

KID=r'results/kidney_bulk/deseq2_out'
MUS=r'results/muscle_bulk/v2/deseq2_out'
CONTRASTS=['age','intv','intv1','intv2','combi_ctrl']  # age vs ctrl, combi vs age, DR vs age, sAct vs age, combi vs ctrl

def load(dirp):
    d={}
    for c in CONTRASTS:
        df=pd.read_csv(f'{dirp}/res_{c}_df.csv')
        d[c]=df
    return d

def build(dirp, sym_from=None):
    d=load(dirp)
    # unify key: gene_id; attach symbol
    base=d['age'][['gene_id']].copy()
    if 'gene_symbol' in d['age'].columns:
        base['sym']=d['age']['gene_symbol']
    elif sym_from is not None:
        base=base.merge(sym_from, on='gene_id', how='left')
    m=base.copy()
    for c in CONTRASTS:
        sub=d[c][['gene_id','log2FoldChange','padj']].rename(columns={'log2FoldChange':f'lfc_{c}','padj':f'padj_{c}'})
        m=m.merge(sub,on='gene_id',how='outer')
    return m, d

def sgn(x):
    return np.sign(x)

def classify(m):
    age_lfc=m['lfc_age']; age_padj=m['padj_age']
    is_deg = age_padj.notna() & (age_padj<0.05) & age_lfc.notna()
    deg=m[is_deg].copy()
    a=sgn(deg['lfc_age'])
    rev_combi = sgn(deg['lfc_intv'])== -a
    rev_dr    = sgn(deg['lfc_intv1'])== -a
    rev_sact  = sgn(deg['lfc_intv2'])== -a
    rev_combi=rev_combi.fillna(False); rev_dr=rev_dr.fillna(False); rev_sact=rev_sact.fillna(False)
    rescued = rev_combi
    norm = deg['lfc_combi_ctrl'].abs()<0.5
    rtn = rescued & norm.fillna(False)
    non_rescued = ~(rev_combi|rev_dr|rev_sact)
    # drivers among RTN
    r=deg[rtn]
    rdr=sgn(r['lfc_intv1'])==-sgn(r['lfc_age']); rdr=rdr.fillna(False)
    rsa=sgn(r['lfc_intv2'])==-sgn(r['lfc_age']); rsa=rsa.fillna(False)
    dr_excl=(rdr&~rsa).sum(); sa_excl=(rsa&~rdr).sum(); dual=(rdr&rsa).sum(); combi_only=(~rdr&~rsa).sum()
    return dict(
        n_genes=len(m), n_deg=int(is_deg.sum()),
        n_rescued=int(rescued.sum()), n_rtn=int(rtn.sum()),
        pct_rescued=round(100*rescued.sum()/is_deg.sum(),1),
        pct_rtn_of_deg=round(100*rtn.sum()/is_deg.sum(),1),
        pct_rtn_of_rescued=round(100*rtn.sum()/rescued.sum(),1),
        n_nonrescued=int(non_rescued.sum()), pct_nonrescued=round(100*non_rescued.sum()/is_deg.sum(),1),
        dr_excl=int(dr_excl), sa_excl=int(sa_excl), dual=int(dual), combi_only=int(combi_only),
        pct_dr=round(100*dr_excl/rtn.sum(),1), pct_sa=round(100*sa_excl/rtn.sum(),1),
        pct_dual=round(100*dual/rtn.sum(),1), pct_combi=round(100*combi_only/rtn.sum(),1),
    ), deg, rtn

def gene_lookup(m, syms):
    out={}
    for s in syms:
        row=m[m['sym']==s] if 'sym' in m.columns else m[m.get('sym','')==s]
        if len(row)==0: out[s]=None; continue
        row=row.iloc[0]
        out[s]=dict(age_lfc=_r(row['lfc_age']),age_padj=_e(row['padj_age']),
                    combi_lfc=_r(row['lfc_intv']),combi_padj=_e(row['padj_intv']),
                    combi_ctrl_lfc=_r(row['lfc_combi_ctrl']),combi_ctrl_padj=_e(row['padj_combi_ctrl']))
    return out
def _r(x):
    return None if pd.isna(x) else round(float(x),3)
def _e(x):
    return None if pd.isna(x) else float(f'{x:.2e}')

def run_once():
    km,_=build(KID)
    sym_map=km[['gene_id','sym']].dropna().drop_duplicates('gene_id')
    mm,_=build(MUS, sym_from=sym_map)
    kc,kdeg,krtn=classify(km)
    mc,mdeg,mrtn=classify(mm)
    # cross-tissue shared RTN by symbol
    ks=set(km.loc[krtn.index[krtn.values] if False else krtn.index,'sym']) if False else None
    krtn_syms=set(kdeg.loc[krtn.index[krtn.values],'sym'].dropna()) if hasattr(krtn,'values') else set()
    # proper: get RTN rows
    def rtn_syms(deg,rtn):
        return deg[rtn]
    krows=kdeg[krtn]; mrows=mdeg[mrtn]
    kset=set(krows['sym'].dropna()); mset=set(mrows['sym'].dropna())
    shared=kset & mset
    # directional consensus + pearson on shared, using aging LFC of both
    kk=krows.set_index('sym'); mmi=mrows.set_index('sym')
    rows=[]
    for g in shared:
        rows.append((g, kk.loc[g,'lfc_age'] if not isinstance(kk.loc[g,'lfc_age'],pd.Series) else kk.loc[g,'lfc_age'].iloc[0],
                        mmi.loc[g,'lfc_age'] if not isinstance(mmi.loc[g,'lfc_age'],pd.Series) else mmi.loc[g,'lfc_age'].iloc[0]))
    sh=pd.DataFrame(rows,columns=['sym','k_age','m_age']).dropna()
    consensus=int((np.sign(sh['k_age'])==np.sign(sh['m_age'])).sum())
    from scipy.stats import pearsonr
    r,p=pearsonr(sh['k_age'],sh['m_age'])
    cross=dict(n_shared=len(shared), pct_of_kidney=round(100*len(shared)/len(kset),1),
               pct_of_muscle=round(100*len(shared)/len(mset),1),
               consensus=f'{consensus}/{len(sh)}', pct_consensus=round(100*consensus/len(sh),1),
               pearson_r=round(float(r),3), pearson_p=float(f'{p:.2e}'))
    # specific genes
    kid_genes=['Tgfb1','Eda2r','Cdkn2a','Cdkn1a','Klotho','Kl','Sirt3','Foxo1','Sod1','Cat','Bax','Lmna',
               'Igf1r','Insr','Braf','Mapk1','Mapk3','Akt1','Akt2','Pten','Prkaa1','Prkab2',
               'Arhgef11','Dbnl','Cfl1','Grap','Arhgap5','Crkl']
    mus_genes=['Tgfb1','Gdf15','Kcnq5','Ptprm','Ryr3','Nav3','Wwox','Raf1','Prkaa2','Bax','Lmna','Ttn','Oxr1','Nrros','Gpx1','Col1a1','Col3a1','Mmp2']
    kg=gene_lookup(km,kid_genes)
    mg=gene_lookup(mm,mus_genes)
    # membership of nutrient/ras genes in RTN sets
    def in_rtn(genes,rows):
        s=set(rows['sym'].dropna()); return {g:(g in s) for g in genes}
    kid_rtn_membership=in_rtn(['Igf1r','Insr','Braf','Mapk1','Mapk3','Akt1','Akt2','Pten','Prkaa1','Prkab2','Tgfb1','Bax','Lmna','Arhgef11','Dbnl','Cfl1','Grap','Arhgap5','Crkl'],krows)
    mus_rtn_membership=in_rtn(['Raf1','Prkaa2','Tgfb1','Bax','Lmna','Kcnq5','Ptprm','Ryr3','Nav3','Wwox','Oxr1','Nrros','Gpx1','Ttn'],mrows)
    return dict(kidney=kc, muscle=mc, cross=cross, kid_genes=kg, mus_genes=mg,
                kid_rtn_membership=kid_rtn_membership, mus_rtn_membership=mus_rtn_membership)

# run 3x consensus
runs=[run_once() for _ in range(3)]
def keyset(d):
    return json.dumps(d, sort_keys=True, default=str)
agree = keyset(runs[0])==keyset(runs[1])==keyset(runs[2])
print('=== 3-RUN CONSENSUS:', 'ALL 3 IDENTICAL' if agree else 'MISMATCH!!!')
print(json.dumps(runs[0], indent=1, default=str))
