# -*- coding: utf-8 -*-
import h5py, numpy as np
def cats(g):
    c=g['categories'][:]; codes=g['codes'][:]
    cc=[x.decode() if isinstance(x,bytes) else x for x in c]
    return np.array([cc[i] if i>=0 else None for i in codes],dtype=object)
h=h5py.File(r'original_data/sn/annotated-aging-soupxed.h5ad','r')
samp=cats(h['obs']['sample'])        # named condition
orig=cats(h['obs']['orig.ident'])    # sg IDs
h.close()
print('KIDNEY atlas crosstab (orig.ident -> sample):')
import collections
m=collections.defaultdict(collections.Counter)
for o,s in zip(orig,samp): m[o][s]+=1
for o in sorted(m):
    top=m[o].most_common(1)[0]
    print(f'  {o} -> {top[0]}   (n={top[1]})')
