#!/usr/bin/env python
# combine for results
from tqdm.auto import tqdm
import os
import pandas as pd 
import sys
import warnings

warnings.filterwarnings("ignore")


from parallel_pandas import ParallelPandas
ParallelPandas.initialize()


pair_fn = sys.argv[1]
db_path = sys.argv[2]
out_fn = sys.argv[3]

pair_df = pd.read_csv(pair_fn, sep='\t')
for (sp, zdf) in tqdm(pair_df.groupby('species')):
    print(sp)
    rmsd = pd.read_hdf(f'{db_path}/{sp}_1.h5', sp).set_index(['gene1', 'gene2'])['lddt'].to_dict()
    zdf['struct_ldo_mdo_lddt'] = zdf.p_apply(lambda x: rmsd.get((x['ldo_gene'], x['mdo_gene']), None), axis=1)

    rmsd = pd.read_hdf(f'{db_path}/{sp}_2.h5', sp).set_index(['gene1', 'gene2'])['lddt'].to_dict()
    zdf['struct_ldo_og1_lddt'] = zdf.p_apply(lambda x: rmsd.get((x['ldo_gene'], x['out_gene']), None) if pd.notna(x['out_gene']) else None, axis=1)

    rmsd = pd.read_hdf(f'{db_path}/{sp}_3.h5', sp).set_index(['gene1', 'gene2'])['lddt'].to_dict()
    zdf['struct_mdo_og1_lddt'] = zdf.p_apply(lambda x: rmsd.get((x['mdo_gene'], x['out_gene']), None) if pd.notna(x['out_gene']) else None, axis=1)

    rmsd = pd.read_hdf(f'{db_path}/{sp}_4.h5', sp).set_index(['gene1', 'gene2'])['lddt'].to_dict()
    zdf['struct_ldo_og2_lddt'] = zdf.p_apply(lambda x: rmsd.get((x['ldo_gene'], x['new_out_gene']), None), axis=1)

    rmsd = pd.read_hdf(f'{db_path}/{sp}_5.h5', sp).set_index(['gene1', 'gene2'])['lddt'].to_dict()
    zdf['struct_mdo_og2_lddt'] = zdf.p_apply(lambda x: rmsd.get((x['mdo_gene'], x['new_out_gene']), None), axis=1)

    rmsd = pd.read_hdf(f'{db_path}/{sp}_6.h5', sp).set_index(['gene1', 'gene2'])['lddt'].to_dict()
    zdf['struct_og1_og2_lddt'] = zdf.p_apply(lambda x: rmsd.get((x['out_gene'], x['new_out_gene']), None), axis=1)

    # need to split the output based on our choice of additional outgroup
    for (ogtype, zzdf) in zdf.groupby('new_out_type'):
        zzdf.to_hdf(out_fn + ogtype + '.h5', sp)
