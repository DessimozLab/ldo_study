#!/usr/bin/env python
'''
    Script to add the expression column to our pairwise tests for a single species
'''
from property_manager import lazy_property
from tqdm.auto import tqdm
import numpy as np
import os
import pandas as pd
import sys


from property_manager import lazy_property
import os



class Expression(object):
    def __init__(self, sp, expression_path='./results/expression/'):
        self.sp = sp
        self._expression_path = expression_path

    @lazy_property
    def genes(self):
        # load the gene order in the expression matrix
        fn = os.path.join(self._expression_path, 'raw', f'{self.sp}_genes.txt')
        if not os.path.isfile(fn):
            raise ValueError(f'No file for species {self.sp}.')
        with open(fn, 'rt') as fp:
            return {j: i for (i, j) in enumerate(map(lambda x: x.rstrip(), fp))}

    @lazy_property
    def raw(self):
        # load the raw expression calls
        fn = os.path.join(self._expression_path, 'raw', f'{self.sp}.npz')
        return np.load(fn)['expr']

    def filtered(self, cutoff):
        # loads the entity-filtered so that we have cutoff % completeness for a particular term
        c = '_'.join(f'{cutoff:.02f}'.split('.'))
        fn = os.path.join(self._expression_path, 'completeness', f'{self.sp}_{c}.npz')
        return np.load(fn)['expr']

    @lazy_property
    def majority_rule(self):
        # load the majority rule expression calls
        fn = os.path.join(self._expression_path, 'majority_rule', f'{self.sp}.npz')
        return np.load(fn)['expr']

    @lazy_property
    def nearest_neighbour(self):
        # load the nearest neighbour expression calls
        fn = os.path.join(self._expression_path, 'nearest_neighbour', f'{self.sp}.npz')
        return np.load(fn)['expr']


def compute_distance(expr, i, j, ignore_missing):
    '''
        Compute hamming distance, making allowance for missing data
    '''
    # load expression
    e_i = expr[i]
    e_j = expr[j]
    n = e_i.shape[0]

    # filter
    z_i = (expr[i] == 0)
    z_j = (expr[j] == 0)

    if not ignore_missing and ((z_i.sum() > 0) or (z_j.sum() > 0)):
        return -1

    elif ignore_missing:
        if (z_i & z_j).sum() > 0:
            # do not accept when both have missing data for a tissue
            return -1

        else:
            f = ~(z_i | z_j)
            e_i = e_i[f]
            e_j = e_j[f]

    if len(e_i) > 0:
        return (e_i != e_j).sum() / n
    else:
        return -1


def compute_expression(df, expr):
    cutoffs = np.arange(0, 1.05, 0.05)
    for cutoff in tqdm(cutoffs):
        c = '_'.join(f'{cutoff:.02f}'.split('.'))

        e = expr.filtered(cutoff)

        df[f'expr_filtered_{c}'] = df.apply(lambda x: compute_distance(expr.filtered(cutoff),
                                                                       x['ldo_g_ii'],
                                                                       x['other_g_ii'],
                                                                       ignore_missing=False),
                                            axis=1)
        df[f'expr_filtered_ignore_{c}'] = df.apply(lambda x: compute_distance(expr.filtered(cutoff),
                                                                              x['ldo_g_ii'],
                                                                              x['other_g_ii'],
                                                                              ignore_missing=True),
                                                   axis=1)

    df['expr_majority_rule'] = df.apply(lambda x: compute_distance(expr.majority_rule,
                                                                   x['ldo_g_ii'],
                                                                   x['other_g_ii'],
                                                                   ignore_missing=False),
                                        axis=1)
    df['expr_nearest_neighbour'] = df.apply(lambda x: compute_distance(expr.nearest_neighbour,
                                                                       x['ldo_g_ii'],
                                                                       x['other_g_ii'],
                                                                       ignore_missing=False),
                                            axis=1)

    return df


def main(args):
    (sp, expression_path, pairs_fn, res_fn) = args[:4]
    # load pairs
    pairs_df = pd.read_csv(pairs_fn, sep='\t')

    # load expr
    expr = Expression(sp, expression_path=expression_path)
    df = compute_expression(pairs_df, expr)
    df.to_csv(res_fn, sep='\t', index=False)


if __name__ == '__main__':
    main(sys.argv[1:])

