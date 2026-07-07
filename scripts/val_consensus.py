# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
KRTN=r'results/Kidney_analysis_results/1_DEGs_tables/kidney_rescued_to_normal.csv'
MRTN=r'results/Muscle_v2_analysis_results/1_DEGs_tables/muscle_v2_rescued_to_normal.csv'
named13=['Crybg3','Eef1a1','Esd','Lgmn','Lztfl1','Man2b1','Minpp1','Ninj1','Nsmce4a','Tars2','Tbc1d1','Ubxn4','Zftraf1']
kr=pd.read_csv(KRTN); mr=pd.read_csv(MRTN)
# age-vs-ctrl LFC column = 'log2FoldChange' in each RTN file (the same LFC the rescue set was built from)
k=kr.dropna(subset=['Symbol']).drop_duplicates('Symbol').set_index('Symbol')['log2FoldChange']
m=mr.dropna(subset=['Symbol']).drop_duplicates('Symbol').set_index('Symbol')['log2FoldChange']
shared=sorted(set(k.index)&set(m.index))
df=pd.DataFrame({'k':k.reindex(shared),'m':m.reindex(shared)})
df['agree']=np.sign(df['k'])==np.sign(df['m'])
print('shared:',len(df),' agree:',int(df['agree'].sum()),' opposite:',int((~df['agree']).sum()))
print('Pearson r:',round(df[['k','m']].corr().iloc[0,1],4))
opp=df[~df['agree']].copy(); opp['min_abs']=opp[['k','m']].abs().min(axis=1)
opp=opp.sort_values('min_abs')
print('\n=== OPPOSITE-DIRECTION genes (my recompute), sorted by closeness to zero ===')
print(opp.round(4).to_string())
my_opp=set(opp.index)
print('\nIn your named-13 but AGREES in my recompute (the culprit):', sorted(set(named13)-my_opp))
print('In my opposite-set but NOT in your named-13:', sorted(my_opp-set(named13)))
