f'''
    Tools to parse the PANTHER trees and estimate relative family rates.
'''
from collections import Counter, defaultdict
from dask.diagnostics import ProgressBar
from enum import Enum
from dendropy import Tree
from functools import lru_cache
from property_manager import cached_property, lazy_property
from scipy.sparse import csr_matrix, linalg
from scipy.spatial.distance import cosine
from tqdm.auto import tqdm
from tqdm.dask import TqdmCallback
import dask.bag as db
import dask.dataframe as dd
import itertools
import numpy as np
import pandas as pd
import re
import sys
import os


class Node(object):
    def __init__(self, n):
        self.n = n

    def get_label(self):
        if self.n.is_leaf():
            return self.n.taxon.label
        else:
            return self.n.label

    def get_parent(self):
        return Node(self.n.parent_node)

    @classmethod
    def find(cls, t, label):
        n = t.find_node_with_label(label)
        if n is None:
            n = t.find_node_with_taxon_label(label)
        return cls(n)


class Taxonomy(object):
    def __init__(self, fn):
        self.t = Tree.get(path=fn, schema='newick')
        for n in self.t.preorder_internal_node_iter():
            n.label = n.annotations['S'].value

        self.labels = {}
        self.i_to_label = []

    @lazy_property
    def terminal_edges(self):
        return {i for (e, i) in self.edges.items() if e[1] in self.taxon_idx}

    @lazy_property
    def taxon_idx(self):
        return set(map(self.get_label_i,
                       map(lambda n: n.taxon.label,
                           self.t.leaf_nodes())))

    @lazy_property
    def valid(self):
        z = set()
        for e in self.t.edges():
            if e.tail_node is not None:
                i = self.get_label_i(e.tail_node.label)
                j = self.get_label_i(e.head_node.taxon.label if e.is_terminal() else e.head_node.label)
                z.add((i, j))
        return z

    @lazy_property
    def edges(self):
        return {e: i for (i, e) in enumerate(self.valid)}

    def get_label_i(self, x):
        if x not in self.labels and x != '':
            self.i_to_label.append(x)
            self.labels[x] = len(self.i_to_label)-1

        return self.labels.get(x, -1)

    def fit_distances(self, dists):
        for e in self.t.edges():
            if e.tail_node is not None:
                i = self.get_label_i(e.tail_node.label)
                j = self.get_label_i(e.head_node.taxon.label if e.is_terminal() else e.head_node.label)
                ij = self.get_path(i, j)
                assert len(ij) == 1
                ij = ij[0]
                d = dists[ij]
                e.length = d

    def get_node(self, i):
        label = self.i_to_label[i]
        return Node.find(self.t, label)

    @lru_cache(None)
    def get_path(self, tail, head):
        '''
            Find path between tail -> head
            NOTE: i is child, j must be on path to root!
        '''
        start = self.i_to_label[head]
        end = self.i_to_label[tail]

        path = []

        n1 = self.get_node(head)
        n2 = n1.get_parent()
        while True:
            path.append(self.edges[(self.labels[n2.get_label()],
                                    self.labels[n1.get_label()])])
            if n2.get_label() == end:
                break

            n1 = n2
            n2 = n1.get_parent()
        return path


class PantherTrees(object):
    #PATTERN = re.compile(r'(?P<an>AN[0-9]*):(?P<sp>[A-Z0-9]*)\|.*UniProtKB=(?P<uniprot>[A-Z0-9]*)')

    def __init__(self, tree_path, taxonomy_fn):
        '''
            Takes as input path to where the PANTHER trees are stored and loads
            them.
            TODO: include expression mapping from bgee. (pass a mapping file for the UniProtKB IDs -> ensembl IDs used in bgee.
        '''
        self.tree_path = tree_path
        self.taxonomy = Taxonomy(taxonomy_fn)

        # solution to the rates / taxonomy branches problem
        self.__solution = None

    @lazy_property
    def _family_files(self):
        z = {}
        for fn in filter(lambda fn: fn.startswith('PTHR') and fn.endswith('.tree'),
                         os.listdir(self.tree_path)):
            fam = int(fn[4:].split('.')[0])
            z[fam] = os.path.join(self.tree_path, fn)
        return z

    @lazy_property
    def uniprot_ids(self):
        return set(map(lambda x: x['UniProtKB'],
                       itertools.chain(*map(lambda x: x.genes.values(),
                                            self.iter_fams()))))

    def get_fam(self, fam_id):
        fn = self._family_files[fam_id]
        with open(fn, 'rt') as fp:
            nwk = fp.readline()
            #genes = {int(m['an'][2:]): {'species': m['sp'],
            #                            'uniprot': m['uniprot']}
            #         for m in map(self.PATTERN.match, fp)}
            genes = {}
            for x in fp:
                x = x.rstrip()
                assert x.startswith('AN') and x[-1] == ';', '{} {}'.format(fam_id, x)
                x = x[:-1] # remove ';'
                
                (an, *y) = x.split(':')
                an_i = int(an[2:])
                y = ':'.join(y).split('|')
                sp = y[0]
                genes[an_i] = {'species': sp}
                for z in y[1:]:
                    (k, *v) = z.split('=')
                    v = '='.join(v)
                    genes[an_i][k] = v                

        return PantherTree(fam_id=fam_id, nwk=nwk, genes=genes)

    def iter_fams(self, progress=True):
        yield from map(self.get_fam,
                       tqdm(self._family_files.keys(),
                            disable=(not progress),
                            desc='For each family'))

    @lazy_property
    def distances(self):
        '''
            Loads distances
        '''
        columns=['fam_id', 'filter_edge',
                 'tail', 'fam_tail_idx', 'fam_tail_ev',
                 'head', 'fam_head_idx', 'fam_head_ev',
                 'length', 'over_duplication']

        def compute_path(x):
            if not x['filter_edge']:
                i = self.taxonomy.get_label_i(x['tail'])
                j = self.taxonomy.get_label_i(x['head'])
                return self.taxonomy.get_path(i, j)

        fam_ids = db.from_sequence(self._family_files)
        print(' - Parsing families', file=sys.stderr)
        with ProgressBar(minimum=1.0, out=sys.stderr):
            df = pd.DataFrame(itertools.chain.from_iterable(fam_ids.map(lambda i:
                                                                        self.get_fam(i).get_distances())),
                                                            columns=columns)
        print(' - Mapping edges', file=sys.stderr)
        df['path'] = df.apply(compute_path, axis=1)
        df['edge'] = df['path'].apply(lambda x: x[0] if (x is not None) and (len(x) == 1) else None)

        return df

    @cached_property
    def matrices(self):
        def family_topology(df):
            indices = []
            indptr = [0]
            for x in df['path'].apply(sorted):
                indices += x
                indptr.append(len(indices))
            # TODO: check dtype
            data = np.ones(len(indices), dtype=np.int8)
            return (csr_matrix((data, indices, indptr)), df.length.values)

        def family_membership(df, fam_ids):
            f2i = {f: i for (i, f) in enumerate(fam_ids)}

            indices = df['fam_id'].apply(f2i.__getitem__)
            indptr = np.arange(0, len(indices)+1, dtype=np.uint32)
            # TODO: check dtype
            data = np.ones(len(indices), dtype=np.int8)
            return (csr_matrix((data, indices, indptr)), f2i)

        def family_rates(df):
            fam_ids = []
            indices = []
            indptr = [0]
            data = []
            obs = []
            for (fam_id, zdf) in df.groupby('fam_id', sort=False):
                (ii, dd) = np.unique(list(itertools.chain.from_iterable(zdf.path)),
                                     return_counts=True)

                indices.append(ii)
                indptr.append(indptr[-1]+len(ii))
                data.append(dd)

                fam_ids.append(fam_id)

                # collect the sum of observations
                obs.append(zdf['length'].sum())

            indices = np.concatenate(indices)
            indptr = np.array(indptr)
            data = np.concatenate(data)
            obs = np.array(obs)

            M = csr_matrix((data, indices, indptr))

            return (M, obs, fam_ids)

        def estimate_tree_dist_vector(df, func):
            df = df[~df['edge'].isna()]
            return getattr(df.groupby('edge')['length'], func)().values

        df = self.distances
        filtered_df = df[~df.filter_edge]
        z = family_rates(filtered_df)

        return {'family_topology': family_topology(filtered_df),             # M, c
                'family_membership': family_membership(filtered_df, z[2]),   # F, 2i
                'family_rates': z,                                           # A, b, fam_ids
                'initial_species': estimate_tree_dist_vector(filtered_df, func='median')  # x0
               }

    def solve_rates(self, start, min_rate=1e-9, eps=1e-9, round_final=3):
        '''
            Identifies relative rates
        '''
        (M, c) = self.matrices['family_topology']
        (A, b, fam_ids) = self.matrices["family_rates"]
        (F, f2i) = self.matrices['family_membership']
        x0 = self.matrices['initial_species']

        b = np.maximum(eps, b)
        c = np.maximum(eps, c)

        x = np.ones_like(x0) if start == 1 else np.maximum(eps, x0)
        r = np.maximum(min_rate, b / (A * x))
        r /= np.mean(r)

        xs = [x]
        rs = [r]
        errs = []
        i = 0
        while True:
            # refine topological distances
            d = c / (F * rs[-1])
            (x1, istop, itn, *_) = linalg.lsqr(M, (c / (F * rs[-1])))
            x1 = np.maximum(eps, x1)
            xs.append(x1)

            # estimate rates
            r = np.maximum(min_rate, b / (A * xs[-1]))
            r /= np.mean(r)
            rs.append(r)
            i += 1

            if len(rs) > 2 and len(xs) > 1:
                dr = cosine(rs[-1], rs[-2])
                dx = cosine(xs[-1], xs[-2])
                l2 = np.linalg.norm((F*r)*(M * xs[-1]) - c, 2)
                l2_normalised = l2 / np.sum(c)
                errs.append((dr, dx, l2_normalised))
                print(i, dr, dx, l2_normalised)
                if dr < eps and dx < eps:
                    break

        if round_final > 0:
            x = np.round(xs[-1], round_final)
            r = np.maximum(min_rate, b / (A * x))
        else:
            x = xs[-1]
            r = rs[-1]

        # cache the result
        self.__solution = {'taxonomy': x, 'rates': r, 'iters': (xs, rs, errs)}
        return self.__solution

    def compute_expected_branch(self, df, start=0):#, eps=1e-3, min_rate=1e-12):
        '''
            Takes a dataframe (formatted as in PantherTrees.distances) and adds
            "expected length"
            NOTE: start=0 indicates using median matching branches to start, start=1 uses all 1.
        '''
        f = df.filter_edge
        z = np.zeros(len(df), dtype=np.float64)
        z[f] = np.nan

        if self.__solution is None:
            ret = self.solve_rates(start=start)#eps=eps, min_rate=min_rate)
        else:
            ret = self.__solution


        f2i = self.matrices['family_membership'][1]

        # Take path & sum edges
        x = df[~f]["path"].apply(lambda x: ret['taxonomy'][x].sum())
        # Load the relevant family rate
        y = df[~f]["fam_id"].apply(lambda x: ret['rates'][f2i[x]])
        z[~f] = (x * y)
        return np.round(z, 3) #(x * y)


class PantherTree(object):
    def __init__(self, fam_id, nwk, genes):
        self.id = fam_id
        self.genes = genes
        self._load(nwk)

    def get_anx(self, n):
        if n.is_leaf():
            i = int(n.taxon.label[2:])
            return (i, self.genes[i]['species'])
        else:
            return (int(n.annotations['ID'].value[2:]), n.label)

    def get_distances(self):
        '''
           Yields distances from family
        '''
        FILTERED_EVENTS = {EvolutionaryEvent.HORIZONTAL_TRANSFER,
                           EvolutionaryEvent.UNKNOWN}
        for e in self.t.preorder_edge_iter():
            if e.tail_node is not None:
                (tail_node, head_node) = (e.tail_node, e.head_node)
                path_length = e.length
                filter_edge_len_2 = (path_length == 2.0)

                if (head_node.event == EvolutionaryEvent.DUPLICATION):
                    # collect paralogous relationships when one level down
                    continue

                over_duplication = (tail_node.event == EvolutionaryEvent.DUPLICATION)
                if over_duplication:
                    # go up until we hit an orthologous relationship
                    while True:
                        # update branch length
                        if tail_node.edge.length is not None:
                            path_length += tail_node.edge.length
                            if not filter_edge_len_2 and tail_node.edge.length == 2.0:
                                filter_edge_len_2 = True
                        tail_node = tail_node.parent_node
                        if (tail_node is None) or (tail_node.event !=
                                                   EvolutionaryEvent.DUPLICATION):
                            break

                    if tail_node is None:
                        # filter out root-level duplications
                        continue

                # filter events / remove tagged poorly-fitting edges
                filter_edge = ((tail_node.event in FILTERED_EVENTS) or
                               (head_node.event in FILTERED_EVENTS) or
                               filter_edge_len_2)

                (an_i, tail) = self.get_anx(tail_node)
                (an_j, head) = self.get_anx(head_node)

                yield (self.id, filter_edge, tail, an_i, tail_node.event.name, head, an_j,
                       head_node.event.name, path_length, over_duplication)

    def get_leafset_taxa(self, i):
        return set(map(lambda n: int(n.taxon.label[2:]),
                       self.nodes[i].leaf_iter()))

    def _load(self, nwk):
        self.t = Tree.get(data=nwk, schema='newick')
        self.t.encode_bipartitions(collapse_unrooted_basal_bifurcation=False)

        # Nodes are labeled "ANx" in preorder in PANTHER trees.
        self.nodes = list(self.t.preorder_node_iter())
        self.n = len(self.t.leaf_nodes())

        self.leaves = set(map(lambda n: int(n.taxon.label[2:]),
                              self.t.leaf_node_iter()))

        for n in self.t.preorder_internal_node_iter():
            n.event = EvolutionaryEvent.identify(n)
            n.label = n.annotations['S'].value

        for n in self.t.leaf_node_iter():
            n.event = EvolutionaryEvent.EXTANT_GENE

        #self._annotate_ldo()

    @lazy_property
    def ldo_groups(self):
        # first, pass up the tree to find distances
        distances = defaultdict(lambda: defaultdict(lambda: np.inf))
        n_genes = Counter()
        for n in self.nodes[:1:-1]:
            label = self.genes[int(n.taxon.label[2:])]['species'] if n.is_leaf() else n.label
            if len(label) > 0:
                distances[n][label] = 0
    
            if n.parent_node.event in {EvolutionaryEvent.SPECIATION, EvolutionaryEvent.DUPLICATION}:
                n_genes[n.parent_node] += (1 if n.is_leaf() else n_genes[n])
                for label in distances[n]:
                    distances[n.parent_node][label] = min(distances[n.parent_node][label], (distances[n][label] + n.edge.length))

        # then pass down the LDO groupings.
        # here we create new groupings when nothing comparable between LHS and RHS. Could instead use gene-weight
        self.t.seed_node.grp = 1
        i = 2

        for n in self.t.preorder_internal_node_iter():
            ch = n.child_nodes()

            if n.event == EvolutionaryEvent.DUPLICATION:
                # LDO == shortest path; others = new group
                # Note: this is performed pairwise and then the shortest absolute path is taken for any shared label between 2 children.
                #       This generalisation is necessary for multifurcating duplication events.
                ch_choice = []
                for (c1, c2) in itertools.combinations(range(len(ch)), 2):
                    # identify shared labels
                    shared_labels = None
                    zd = defaultdict(lambda: np.inf)
                    zn = {}
                    #for (k, c) in enumerate(ch):
                    for (k, c) in [(c1, ch[c1]), (c2, ch[c2])]:
                        if c.event in {EvolutionaryEvent.SPECIATION, EvolutionaryEvent.DUPLICATION, EvolutionaryEvent.EXTANT_GENE}:
                            for (d, e) in distances[c].items():
                                if e + c.edge.length < zd[d]:
                                    zd[d] = e + c.edge.length
                                    zn[d] = k

                            if shared_labels is None:
                                shared_labels = set(distances[c].keys())
                            else:
                                shared_labels &= set(distances[c].keys())

                    if shared_labels is not None and len(shared_labels) > 0:
                        # identify shortest label
                        min_label = None
                        for x in shared_labels:
                            if min_label is None or zd[x] < zd[min_label]:
                                min_label = x
                        ch_choice.append((zn[min_label], zd[min_label]))
                        #j = zn[min_label]
    
                    #else:
                    #    # len(shared_labels) == 0:
                    #    j = None
                    #    '''# ALTERNATIVE: take heaviest side
                    #    for (k, c) in enumerate(ch):
                    #        if c.event in {EvolutionaryEvent.SPECIATION, EvolutionaryEvent.DUPLICATION, EvolutionaryEvent.EXTANT_GENE}:
                    #            if j is None or n_genes[c] > n_genes[ch[j]]:
                    #                j = k
                    #    '''
                j = min(ch_choice, key=lambda x: x[1])[0] if len(ch_choice) > 0 else None
                for (k, c) in enumerate(ch):
                    if j != k:
                        c.grp = i
                        i += 1
                    else:
                        # chosen!
                        c.grp = n.grp

            elif n.event in {EvolutionaryEvent.SPECIATION, EvolutionaryEvent.HORIZONTAL_TRANSFER}:
                # pass down set.
                for c in ch:
                    if c.event in {EvolutionaryEvent.SPECIATION, EvolutionaryEvent.DUPLICATION, EvolutionaryEvent.EXTANT_GENE}:
                        c.grp = n.grp
                    else:
                        c.grp = i
                        i += 1

            elif n.event == EvolutionaryEvent.UNKNOWN:
                # set all to new group
                for c in ch:
                    c.grp = i
                    i += 1

        ldo_groups = defaultdict(dict)
        for n in self.t.leaf_nodes():
            an_i = int(n.taxon.label[2:])
            # collect the group -> {sp: uniprot}
            ldo_groups[n.grp][self.genes[an_i]['species']] = self.genes[an_i]['UniProtKB']

        return ldo_groups


class UnknownEvolutionaryEvent(Exception):
    pass

class EvolutionaryEvent(Enum):
    UNKNOWN = -1
    EXTANT_GENE = 0
    SPECIATION = 1
    DUPLICATION = 2
    HORIZONTAL_TRANSFER = 3

    @classmethod
    def identify(cls, n):
        if n.is_leaf():
            return cls.EXTANT_GENE
        elif cls.is_orthologGroup_node(n):
            return cls.SPECIATION
        elif cls.is_paralogGroup_node(n):
            return cls.DUPLICATION
        elif cls.is_hgt_node(n):
            return cls.HORIZONTAL_TRANSFER
        elif cls.is_unknown_evolutionary_node(n):
            return cls.UNKNOWN
        else:
            raise UnknownEvolutionaryEvent('Encountered unknown event type in node:',
                                           Node(n).get_label(),
                                           n.annotations['Ev'].value)

    @staticmethod
    def is_unknown_evolutionary_node(n):
        return (n.annotations['Ev'].value == 'UNK')

    @staticmethod
    def is_orthologGroup_node(n):
        return (n.annotations['Ev'].value == '0>1')

    @staticmethod
    def is_paralogGroup_node(n):
        return (n.annotations['Ev'].value == '1>0')

    @staticmethod
    def is_hgt_node(n):
        return (n.annotations['Ev'].value == '0>0')
