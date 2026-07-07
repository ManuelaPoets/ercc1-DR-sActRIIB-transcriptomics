# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
KBASE=r'results/Kidney_analysis_results/5_Cellchat_LR_analysis'
MBASE=r'results/Muscle_v2_analysis_results/5_Cellchat_LR_analysis'
CONDS=['ctrl','age','DR','sAct','combi']

def load(base, suffix=''):
    d={}
    for c in CONDS:
        d[c]=pd.read_csv(f'{base}/{c}_lr_interactions_reversal_genes{suffix}.csv')
    return d

def report(name, d, pathways):
    print(f'===== {name} =====')
    print('  total interaction prob:', {c:round(float(d[c]["prob"].sum()),1) for c in CONDS})
    print('  interaction count     :', {c:len(d[c]) for c in CONDS})
    pc='pathway_name'
    for pw in pathways:
        vals={c:round(float(d[c].loc[d[c][pc]==pw,'prob'].sum()),3) for c in CONDS}
        # % change ctrl->age
        ch = (f'  ({100*(vals["age"]-vals["ctrl"])/vals["ctrl"]:+.0f}% ctrl→age)' if vals['ctrl'] else '')
        print(f'    {pw:9}: {vals}{ch}')

# 3-run consensus (deterministic)
for run in range(3):
    kd=load(KBASE); md=load(MBASE,'_v2')
    if run==0:
        report('KIDNEY (3.5.1)', kd, ['VEGF','EGF','SEMA3','EPHA','EPHB','PDGF','MIF','NEGR'])
        report('MUSCLE v2 (3.5.2)', md, ['PTPRM','APP','PROS','TGFb','NECTIN'])
import json
def snap():
    kd=load(KBASE); md=load(MBASE,'_v2')
    return json.dumps({c:[round(float(kd[c]['prob'].sum()),3),len(kd[c]),round(float(md[c]['prob'].sum()),3),len(md[c])] for c in CONDS})
print('\n3-RUN CONSENSUS:', 'IDENTICAL' if snap()==snap()==snap() else 'MISMATCH')
