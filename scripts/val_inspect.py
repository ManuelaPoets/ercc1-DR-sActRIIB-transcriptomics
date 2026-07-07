# -*- coding: utf-8 -*-
import h5py, sys
def cats(o):
    try:
        c=o['categories'][:]
        return [x.decode() if isinstance(x,bytes) else x for x in c]
    except Exception as e:
        return None
for f,name in [(r'results/sn/annotated-aging-soupxed.h5ad','KIDNEY'),
               (r'results/sn/annotated-muscle-soupxed.h5ad','MUSCLE')]:
    print('===================='*1, name, f)
    h=h5py.File(f,'r')
    print('top keys:', list(h.keys()))
    print('obs keys:', list(h['obs'].keys()))
    print('X type:', type(h['X']).__name__, (list(h['X'].keys()) if isinstance(h['X'],h5py.Group) else h['X'].shape))
    if 'raw' in h:
        print('raw keys:', list(h['raw'].keys()))
        if 'X' in h['raw']:
            rx=h['raw']['X']
            print('  raw/X:', type(rx).__name__, (list(rx.keys()) if isinstance(rx,h5py.Group) else rx.shape))
        if 'var' in h['raw']:
            print('  raw/var keys:', list(h['raw']['var'].keys()))
    if 'layers' in h: print('layers:', list(h['layers'].keys()))
    print('var keys:', list(h['var'].keys()))
    # var index / gene symbols length
    for vk in ['_index','features','gene_symbols','gene_ids']:
        if vk in h['var']:
            print('  var[%s] n=%d sample=%s'%(vk, h['var'][vk].shape[0], [x.decode() if isinstance(x,bytes) else x for x in h['var'][vk][:3]]))
    if 'raw' in h and 'var' in h['raw']:
        for vk in ['_index','features','gene_symbols','gene_ids']:
            if vk in h['raw']['var']:
                print('  raw/var[%s] n=%d'%(vk, h['raw']['var'][vk].shape[0]))
    for k in h['obs'].keys():
        o=h['obs'][k]
        if isinstance(o,h5py.Group) and 'categories' in o:
            cc=cats(o)
            if cc is not None:
                if len(cc)<=40: print('  obs[%s] (%d):'%(k,len(cc)), cc)
                else: print('  obs[%s]: %d cats'%(k,len(cc)))
    h.close()
