# -*- coding: utf-8 -*-
"""Regenerate Supplementary Table S3.8 — cross-tissue common rescued-to-normal genes.
Built from the canonical v2 RTN tables, intersected by gene symbol."""
import pandas as pd, numpy as np
from scipy.stats import pearsonr

KRTN=r'results/Kidney_analysis_results/1_DEGs_tables/kidney_rescued_to_normal.csv'
MRTN=r'results/Muscle_v2_analysis_results/1_DEGs_tables/muscle_v2_rescued_to_normal.csv'
OUT=r'results/common_rescued_genes_kidney_muscle.csv'

k=pd.read_csv(KRTN).rename(columns={'Unnamed: 0':'gene_id'})
m=pd.read_csv(MRTN).rename(columns={'Unnamed: 0':'gene_id'})
k=k.dropna(subset=['Symbol']).drop_duplicates('Symbol')
m=m.dropna(subset=['Symbol']).drop_duplicates('Symbol')

shared=sorted(set(k['Symbol']) & set(m['Symbol']))
ki=k.set_index('Symbol'); mi=m.set_index('Symbol')

rows=[]
for g in shared:
    kr=ki.loc[g]; mr=mi.loc[g]
    k_id=kr['gene_id']; m_id=mr['gene_id']
    kl=float(kr['log2FoldChange']); ml=float(mr['log2FoldChange'])
    concord = np.sign(kl)==np.sign(ml)
    rows.append({
        'Symbol': g,
        'Ensembl_ID': k_id,
        'Ensembl_ID_muscle': (m_id if m_id!=k_id else ''),  # filled only if it differs (dual-ID genes)
        'Kidney_log2FC_aging': round(kl,3),
        'Kidney_padj_aging': float(f"{kr['padj']:.2e}"),
        'Muscle_log2FC_aging': round(ml,3),
        'Muscle_padj_aging': float(f"{mr['padj']:.2e}"),
        'Direction': 'Concordant' if concord else 'Discordant',
        'Kidney_log2FC_combi_vs_ctrl': round(float(kr['log2FC_combi_ctrl']),3),
        'Muscle_log2FC_combi_vs_ctrl': round(float(mr['log2FC_combi_ctrl']),3),
    })
df=pd.DataFrame(rows)
# sort: concordant first, then by |kidney aging LFC| descending
df['_k']=df['Direction'].map({'Concordant':0,'Discordant':1})
df['_a']=df['Kidney_log2FC_aging'].abs()
df=df.sort_values(['_k','_a'],ascending=[True,False]).drop(columns=['_k','_a']).reset_index(drop=True)
df.insert(0,'Rank',range(1,len(df)+1))
df.to_csv(OUT,index=False)

n=len(df); conc=(df['Direction']=='Concordant').sum(); disc=n-conc
r,p=pearsonr(df['Kidney_log2FC_aging'],df['Muscle_log2FC_aging'])
dual=df[df['Ensembl_ID_muscle']!='']['Symbol'].tolist()
print(f'S3.8 written: {OUT}')
print(f'n shared (by symbol) = {n}')
print(f'concordant = {conc}/{n} ({100*conc/n:.1f}%) | discordant = {disc}')
print(f'Pearson r (aging LFC, kidney vs muscle) = {r:.3f} (p = {p:.2e})')
print(f'dual-Ensembl-ID genes = {dual}')
disc_genes=sorted(df.loc[df['Direction']=='Discordant','Symbol'])
print('discordant genes =', disc_genes)
print()
print('=== full table ===')
with pd.option_context('display.max_rows',None,'display.width',200):
    print(df.to_string(index=False))
