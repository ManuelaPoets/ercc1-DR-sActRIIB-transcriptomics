# -*- coding: utf-8 -*-
import pandas as pd, numpy as np, json
from scipy.stats import pearsonr
KRTN=r'results/Kidney_analysis_results/1_DEGs_tables/kidney_rescued_to_normal.csv'
MRTN=r'results/Muscle_v2_analysis_results/1_DEGs_tables/muscle_v2_rescued_to_normal.csv'
KIMP=r'results/Kidney_analysis_results/3_Effect_comparison_tables/kidney_intervention_impact_comparison.csv'
MIMP=r'results/Muscle_v2_analysis_results/3_Effect_comparison_tables/muscle_v2_intervention_impact_comparison.csv'

def run():
    out={}
    kr=pd.read_csv(KRTN); mr=pd.read_csv(MRTN)
    # 1. threshold check
    out['kidney_RTN_rows']=len(kr); out['muscle_RTN_rows']=len(mr)
    out['kidney_max_abs_combi_ctrl']=round(float(kr['log2FC_combi_ctrl'].abs().max()),4)
    out['muscle_max_abs_combi_ctrl']=round(float(mr['log2FC_combi_ctrl'].abs().max()),4)
    out['kidney_n_with_absCC_0.4to0.5']=int(((kr['log2FC_combi_ctrl'].abs()>=0.4)&(kr['log2FC_combi_ctrl'].abs()<0.5)).sum())
    # 2 & 3. drivers + restoration from impact tables
    for name,imp in [('kidney',pd.read_csv(KIMP)),('muscle',pd.read_csv(MIMP))]:
        dr=imp['Intervention_DR_Effect']=='Rescue'; sa=imp['Intervention_sAct_Effect']=='Rescue'
        d=dict(n=len(imp),dr_excl=int((dr&~sa).sum()),sa_excl=int((sa&~dr).sum()),
               dual=int((dr&sa).sum()),combi_only=int((~dr&~sa).sum()))
        d.update({'pct_dr':round(100*d['dr_excl']/d['n'],1),'pct_sa':round(100*d['sa_excl']/d['n'],1),
                  'pct_dual':round(100*d['dual']/d['n'],1),'pct_combi':round(100*d['combi_only']/d['n'],1)})
        out[name+'_drivers']=d
        dual=imp[dr&sa]
        out[name+'_restore_dual']=dict(
            DR_med=round(float(dual['Restoration_DR'].median()),3),
            sAct_med=round(float(dual['Restoration_sAct'].median()),3),
            pct_DR_gt_sAct=round(100*float((dual['Restoration_DR']>dual['Restoration_sAct']).mean()),1),
            n_dual=len(dual))
        # also medians over ALL rtn
        out[name+'_restore_all']=dict(DR_med=round(float(imp['Restoration_DR'].median()),3),
            sAct_med=round(float(imp['Restoration_sAct'].median()),3))
    # 4. cross-tissue shared by Symbol
    ks=kr.dropna(subset=['Symbol']); ms=mr.dropna(subset=['Symbol'])
    kset=set(ks['Symbol']); mset=set(ms['Symbol']); shared=kset&mset
    # also by Ensembl id (first col)
    kid_ids=set(kr['Unnamed: 0']); mid=set(mr['Unnamed: 0']); shared_id=kid_ids&mid
    kage=ks.drop_duplicates('Symbol').set_index('Symbol')['log2FoldChange']
    mage=ms.drop_duplicates('Symbol').set_index('Symbol')['log2FoldChange']
    rows=[(g,kage[g],mage[g]) for g in shared if g in kage.index and g in mage.index]
    sh=pd.DataFrame(rows,columns=['g','k','m']).dropna()
    cons=int((np.sign(sh['k'])==np.sign(sh['m'])).sum()); r,p=pearsonr(sh['k'],sh['m'])
    out['cross']=dict(shared_by_symbol=len(shared),shared_by_ensembl=len(shared_id),
        kidney_valid_sym=len(kset),muscle_valid_sym=len(mset),
        pct_of_kidney=round(100*len(shared)/len(kset),1),pct_of_muscle=round(100*len(shared)/len(mset),1),
        consensus=f'{cons}/{len(sh)}',pct_consensus=round(100*cons/len(sh),1),
        pearson_r=round(float(r),3),pearson_p=float(f'{p:.2e}'))
    return out

runs=[json.dumps(run(),sort_keys=True) for _ in range(3)]
print('3-RUN CONSENSUS:', 'IDENTICAL' if runs[0]==runs[1]==runs[2] else 'MISMATCH')
print(json.dumps(json.loads(runs[0]),indent=1))
