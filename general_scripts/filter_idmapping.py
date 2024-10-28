#!/usr/bin/env python
'''
    Filter the UniProt idmapping.dat.gz file for particular UniProt entries / xref IDs.
    -- Alex Warwick Vesztrocy, March 2022
'''
from tqdm.auto import tqdm
import pandas as pd
import sys


def read_filter(fn):
    with open(fn, 'rt') as fp:
        return set(map(lambda x: x.rstrip(), fp.readlines()))


#def filter_idmapping(fn, uniprot_ids, xref_ids, out_fp):
def filter_idmapping(fn, xref_ids, out_fp):
    for df in tqdm(pd.read_csv(fn, sep='\t', names=['uniprot_id', 'xref_type', 'xref_id'], chunksize=1e7)):
        #df = df[df['uniprot_id'].isin(uniprot_ids)]
        f1 = df['xref_id'].isin(xref_ids)
        f2 = df['xref_id'].apply(lambda x: (x not in xref_ids) and not pd.isna(x) and (x.split('.')[0] in xref_ids))
        df.loc[f2, 'xref_id'] = df.loc[f2,'xref_id'].apply(lambda x: x.split('.')[0])
        df = df[f1 | f2]

        df.to_csv(out_fp, sep='\t', index=False, header=False)


def main(args):
    #(idmap_fn, uniprot_fn, xref_fn) = args[:3]
    (idmap_fn, xref_fn) = args[:2]

    #uniprot_ids = read_filter(uniprot_fn)
    xref_ids = read_filter(xref_fn)
    #filter_idmapping(idmap_fn, uniprot_ids, xref_ids, out_fp=sys.stdout)
    filter_idmapping(idmap_fn, xref_ids, out_fp=sys.stdout)


if __name__ == '__main__':
    main(sys.argv[1:])
