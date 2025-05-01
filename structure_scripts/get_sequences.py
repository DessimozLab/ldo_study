#!/usr/bin/env python
from tqdm.auto import tqdm
import json
import numpy as np
import os
import pandas as pd
import urllib3
import sys
import tables


def get_url(uniprot_id):
    return f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"


class SpeciesDatabase(object):
    # this would be more efficient to store as a buffer obj with a single table
    def __init__(self, sp, db_path):
        self.db = tables.open_file(os.path.join(db_path, sp+'.h5'), 'w')

    def add(self, id, seq):
        self.db.create_carray('/', id, tables.Atom.from_dtype(seq.dtype), obj=seq, filters=tables.Filters(complevel=6, complib='blosc'))

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

    response = http.request("GET", get_url(uniprot_id), headers={'accept': 'application/json'})
    if response.status == 200:
        x = json.loads(response.data)
        if 'sequence' in x:
            seq = np.frombuffer(x['sequence']['value'].encode('ascii'), dtype='S1')
            db.add(uniprot_id, seq)
            continue

    failed.append(uniprot_id)

db.db.close()


print(' - ', sp, len(failed), 'failed')
for x in failed:
    print('   - ', x)
