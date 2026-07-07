# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
# canonical res_age WITH symbols (what Results used), both tissues
KA=r'results/Kidney_analysis_results/1_DEGs_tables/kidney_res_age_results_with_symbols.csv'
MA=r'results/Muscle_v2_analysis_results/1_DEGs_tables/muscle_v2_res_age_results_with_symbols.csv'
# also combi_ctrl for rescue check
KCC=r'original_data/bulk/kidney_bulk/deseq2_out/res_combi_ctrl_df.csv'
MCC=r'original_data/bulk/muscle_bulk/v2/deseq2_out/res_combi_ctrl_df.csv'
ka=pd.read_csv(KA); ma=pd.read_csv(MA)
print('kidney with_symbols cols:',list(ka.columns)[:8])
print('muscle with_symbols cols:',list(ma.columns)[:8])
def symcol(df):
    for c in ['Symbol','gene_symbol','symbol']:
        if c in df.columns: return c
    return None
def lfccol(df):
    for c in ['log2FoldChange','log2FoldChange_age','lfc']:
        if c in df.columns: return c
    return None
ks,kl=symcol(ka),lfccol(ka); ms,ml=symcol(ma),lfccol(ma)
def get(df,sc,lc,sym):
    r=df[df[sc]==sym]
    if not len(r): return None
    return round(float(r.iloc[0][lc]),3), float(f"{r.iloc[0]['padj']:.2e}")
genes=['Eda2r','Gdf15','Cdkn2a','Cdkn1a','Kl','Klotho','Sirt3','Foxo1','Sod1','Cat','Nfe2l2']
print('%-10s %-30s %-30s'%('gene','KIDNEY age (lfc,padj)','MUSCLE age (lfc,padj)'))
for g in genes:
    k=get(ka,ks,kl,g); m=get(ma,ms,ml,g)
    print(f'{g:10} {str(k):30} {str(m):30}')
