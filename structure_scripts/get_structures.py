#!/usr/bin/env python
from tqdm.auto import tqdm
import numpy as np
import os
import pandas as pd
import urllib3
import sys
import tables


def get_url(uniprot_id):
    return f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'


class SpeciesDatabase(object):
    # this would be more efficient to store as a buffer obj with a single table
    def __init__(self, sp, db_path):
        self.db = tables.open_file(os.path.join(db_path, sp+'.h5'), 'w')

    def add(self, id, pdb):
        self.db.create_carray('/', id, tables.Atom.from_dtype(pdb.dtype), obj=pdb, filters=tables.Filters(complevel=6, complib='blosc'))

    def check(self, id):
        return (id in self.db.root)


sp_fn = sys.argv[1]
sp = sys.argv[2]
db_path = sys.argv[3]

uniprot_ids = list(pd.read_csv(sp_fn, sep='\t', names=['uniprot_id'])['uniprot_id'])

db = SpeciesDatabase(sp, db_path)

failed = []
http = urllib3.PoolManager()
for uniprot_id in tqdm(uniprot_ids):
    if db.check(uniprot_id):
        # already in db
        continue

    response = http.request("GET", get_url(uniprot_id))
    if response.status == 200:
        db.add(uniprot_id, np.frombuffer(response.data, dtype='S1'))
    else:
        failed.append(uniprot_id)

db.db.close()


print(' - ', sp, len(failed), 'failed')
for x in failed:
    print('   - ', x)
