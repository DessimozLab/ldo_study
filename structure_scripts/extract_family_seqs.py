#!/usr/bin/env python
# script to extract sequences for a particular panther family
import pandas as pd
import os
import subprocess
import sys
import tables
import tempfile
import parasail


class SeqDB(object):
    def __init__(self, path):
        self.path = path
        self.dbs = {}

    def get_seq(self, sp, uniprot_id):
        sp_db = self.get_db(sp)
        if uniprot_id in sp_db.root:
            return sp_db.root[uniprot_id][:].tobytes().decode('ascii')

    def get_db(self, sp):
        if sp not in self.dbs:
            self.dbs[sp] = tables.open_file(os.path.join(self.path, f'{sp}.h5'), 'r')
        return self.dbs[sp]

    def close(self):
        for db in self.dbs.values():
            db.close()


df_fn = sys.argv[1]
fam_id = int(sys.argv[2])
seq_path = sys.argv[3]

df = pd.read_csv(df_fn, sep='\t')
df = df[df.fam == int(fam_id)]

seq_db = SeqDB(seq_path)
for (_, r) in df.iterrows():
    s = seq_db.get_seq(r['species'], r['gene_id'])
    if s is not None:
        print('>{}'.format(r['gene_id']))
        while len(s):
            print(s[:60])
            s = s[60:]

seq_db.close()
