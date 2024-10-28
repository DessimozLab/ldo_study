#!/usr/bin/env python
# align for each pair which we have the struct pred for using USalign
from functools import lru_cache
from tqdm.auto import tqdm
import pandas as pd
import os
import subprocess
import sys
import tables
import tempfile


header = ['fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','lddt']


class StructDB(object):
    def __init__(self, path):
        self.path = path
        self.dbs = {}

    def get_struct(self, sp, uniprot_id):
        sp_db = self.get_db(sp)
        if uniprot_id in sp_db.root:
            return sp_db.root[uniprot_id][:].tobytes()

    def get_db(self, sp):
        if sp not in self.dbs:
            self.dbs[sp] = tables.open_file(os.path.join(self.path, f'{sp}.h5'), 'r')
        return self.dbs[sp]

    def close(self):
        for db in self.dbs.values():
            db.close()


def align(struct_db, sp1, uniprot1, sp2, uniprot2):
    pdb1 = struct_db.get_struct(sp1, uniprot1)
    pdb2 = struct_db.get_struct(sp2, uniprot2)
    if pdb1 is not None and pdb2 is not None:
        with tempfile.NamedTemporaryFile() as pdb1_fp, \
                tempfile.NamedTemporaryFile() as pdb2_fp, \
                tempfile.NamedTemporaryFile() as tmp_result_fp, \
                tempfile.TemporaryDirectory() as tmp_path:
            pdb1_fp.write(pdb1)
            pdb1_fp.flush()
            pdb2_fp.write(pdb2)
            pdb2_fp.flush()

            res = subprocess.run(['foldseek', 'easy-search', pdb1_fp.name, pdb2_fp.name, tmp_result_fp.name, tmp_path, '--threads', '1', '--exhaustive-search', '--format-output',  ','.join(header)], capture_output=True)
            pbar.update()
            with open(tmp_result_fp.name, 'rt') as fp:
                x = fp.readline().rstrip().split('\t')
            return [uniprot1, uniprot2] + x

    # default return
    pbar.update()
    return [uniprot1, uniprot2] + [''] * len(header)

pair_fn = sys.argv[1]
sp = sys.argv[2]
struct_path = sys.argv[3]
out_fn = sys.argv[4]
part = int(sys.argv[5])

pair_df = pd.read_csv(pair_fn, sep='\t')
pair_df = pair_df[pair_df.species == sp]
struct_db = StructDB(struct_path)

has_out_f = pair_df.out_gene.notna()

if part == 1:
    pbar = tqdm(total=len(pair_df))
    x = list(pair_df.apply(lambda x: align(struct_db, x['species'], x['ldo_gene'], x['species'], x['mdo_gene']), axis=1))
elif part == 2:
    pbar = tqdm(total=has_out_f.sum())
    x = list(pair_df[has_out_f].apply(lambda x: align(struct_db, x['species'], x['ldo_gene'], x['outgroup_species'], x['out_gene']), axis=1))
elif part == 3:
    pbar = tqdm(total=has_out_f.sum())
    x = list(pair_df[has_out_f].apply(lambda x: align(struct_db, x['species'], x['mdo_gene'], x['outgroup_species'], x['out_gene']), axis=1))
else:
    raise ValueError('unknown part {}'.format(part))

res_df = pd.DataFrame(x, columns=['gene1', 'gene2'] + header)

pbar.close()
struct_db.close()
res_df.to_hdf(out_fn, sp, complevel=9, complib='zlib')
