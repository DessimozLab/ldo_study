#!/usr/bin/env python
'''
    Convert h5 results file to TSV
'''
import pandas as pd
import sys

db = pd.HDFStore(sys.argv[1], 'r')
for sp in map(lambda x: x._v_name, db.root):
    df = pd.read_hdf(db, sp)
    df.to_csv(sys.stdout, sep='\t', index=False)

db.close()
