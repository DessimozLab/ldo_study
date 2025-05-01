#!/usr/bin/env python
'''
Script to generate extra outgroup (og2) that is roughly same distance for LDO-og2 and MDO-og1.
Options:
- 1 use tree distance
- 2 use percentage sequence identity

Note: og1 was chosen as the closest species outgroup, so is not necessarily the closest in the gene tree.
'''
from tqdm.auto import tqdm
import gzip
import pandas as pd
import sys

from lib.PantherParser import PantherTrees, EvolutionaryEvent


def find_nodes(fam, **ii):
    nodes = list(fam.t.preorder_node_iter())
    res = {}
    for (k,i) in ii.items():
        n = nodes[i]
        assert fam.get_anx(n)[0] == i
        res[k] = n
    return res


def compute_distance(fam, M, g1, g2):
    x = find_nodes(fam, g1=g1, g2=g2)
    return M.distance(x['g1'].taxon, x['g2'].taxon)
    

def compute_candidates(align_res_path, fam_id, df):
    fam = pt.get_fam(fam_id)
    uniprot2anx = {v['UniProtKB']: k for (k, v) in fam.genes.items()}
    g2grp = {g: grp for (grp, v) in fam.ldo_groups.items() for (sp, g) in v.items()}
    M = fam.t.phylogenetic_distance_matrix()

    # load pident
    pident_df = pd.read_csv(os.path.join(align_res_path, '{}.tsv.gz'.format(fam.id)), sep='\t')
    not_obsolete = set(pident_df.g1) | set(pident_df.g2)
    swap = pident_df['g1'] > pident_df['g2']
    pident_df = pd.concat((pident_df[~swap], pident_df[swap].rename(columns={'g1': 'g2', 'g2': 'g1'}))).set_index(['g1', 'g2'])
    
    def get_pident(g1, g2):
        return (pident_df.loc[g1,g2] if g1 <= g2 else pident_df.loc[g2,g1])['pident']
        
    for (_, r) in df.iterrows():
        ldo_gene = r['ldo_gene']
        mdo_gene = r['mdo_gene']
        out_gene = r['out_gene']
        out_sp = r['outgroup_species']
        x = find_nodes(fam, ldo=r.ldo_fam_head_idx, mdo=r.mdo_fam_head_idx)

        if ldo_gene not in not_obsolete:
            continue
        if mdo_gene not in not_obsolete:
            continue
        if out_gene not in not_obsolete:
            continue
        if x['ldo'].parent_node != x['mdo'].parent_node:
            # loss events
            continue

        grp = g2grp[ldo_gene]
        if out_gene not in g2grp or grp != g2grp[out_gene]:
            # outgroup situations that we should filter out
            continue

        dup_node = x['ldo'].parent_node
        assert (EvolutionaryEvent.identify(dup_node) == EvolutionaryEvent.DUPLICATION)
        
        # find a further out gene that is in the same group as our LDO / outgroup and similar distance as outgroup is from MDO.
        in_species = set(map(lambda x: x[1], map(fam.get_anx, dup_node.leaf_nodes())))
        sp_choices = set(fam.ldo_groups[grp].keys()) - {out_sp} - in_species
        if len(sp_choices) == 0:
            continue
        
        # MDO - out choice
        mdo_out_dist = compute_distance(fam, M, uniprot2anx[mdo_gene], uniprot2anx[out_gene])
        mdo_out_pident = get_pident(mdo_gene, out_gene)

        ldo_mdo_pident = get_pident(ldo_gene, mdo_gene)
       
        # identify possible choices 
        choices = pd.DataFrame(map(lambda x: [x[0], x[1], compute_distance(fam, M, uniprot2anx[ldo_gene], uniprot2anx[x[1]])],
                                   map(lambda sp: (sp, fam.ldo_groups[grp][sp]), 
                                       sp_choices)),
                               columns=['species', 'new_out_gene', 'dist_to_ldo'])
        choices = choices[choices.new_out_gene.isin(not_obsolete)]
        if len(choices) == 0:
            continue
            
        choices['dist_diff'] = choices['dist_to_ldo'] - mdo_out_dist
        choices['abs_dist_diff'] = choices['dist_diff'].abs()
        choices['pident_to_ldo'] = choices['new_out_gene'].apply(lambda x: get_pident(ldo_gene, x))
        choices['pident_diff'] = choices['pident_to_ldo'] - mdo_out_pident
        choices['abs_pident_diff'] = choices['pident_diff'].abs()

        # choose min abs dist diff
        r = list(r)
        h = ['species', 'new_out_gene', 'dist_to_ldo', 'dist_diff', 'pident_to_ldo', 'pident_diff']
        
        dist_idx = choices['abs_dist_diff'].argmin()
        dist_choice = choices.iloc[dist_idx]
        yield r + ['min_abs_dist'] + list(dist_choice[h]) + [mdo_out_dist, mdo_out_pident, ldo_mdo_pident]
        
        # choose min abs pident diff
        pident_idx = choices['abs_pident_diff'].argmin()
        pident_choice = choices.iloc[pident_idx]
        yield r + ['min_abs_pident'] + list(pident_choice[h]) + [mdo_out_dist, mdo_out_pident, ldo_mdo_pident]

# runner
res_path = sys.argv[1]
align_res_path = sys.argv[2]
tests_fn = sys.argv[3]
my_job_id = int(sys.argv[4])
total_jobs = int(sys.argv[5])

tests_df = pd.read_csv(tests_fn, sep='\t')

h = ['species', 'new_out_gene', 'dist_to_ldo', 'dist_diff', 'pident_to_ldo', 'pident_diff', 'mdo_out_dist', 'mdo_out_pident', 'ldo_mdo_pident']

os.makedirs(res_path, exist_ok=True)

with gzip.open(os.path.join(res_path, 'res{:04d}_{:04d}.tsv.gz'.format(my_job_id, total_jobs)), 'wt') as fp:
    i = 0
    for (fam_id, zdf) in tests_df.groupby('fam_id', sort=False):
        if i%total_jobs == my_job_id:
            z = pd.DataFrame(compute_candidates(fam_id, zdf), columns=list(tests_df.keys()) + ['new_out_type'] + h)
            z.to_csv(fp, sep='\t', index=False)
            fp.flush()
        i += 1
