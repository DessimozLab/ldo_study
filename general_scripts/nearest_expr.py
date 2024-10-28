#!/usr/bin/env python
'''
    Identify closest genes (using known expression profile) and set unknown as most common amongst
    these.
'''
#from sklearn.experimental import enable_iterative_imputer
#from sklearn.impute import IterativeImputer
from fancyimpute import BiScaler
from tqdm.auto import tqdm
import numpy as np
import sys


def nearest_neighbours(expr, cutoff=0.95):
    x = expr.astype(np.float64)
    x[x == 0] = np.nan

    y = BiScaler(min_value=-1, max_value=1).fit_transform(x)
    #y = IterativeImputer(missing_values=0,
    #                     initial_strategy='most_frequent',
    #                     min_value=-1,
    #                     max_value=1,
    #                     tol=1
    #                    ).fit_transform(expr)
    y[y<0] = -1
    y[y>0] = 1
    y = y.astype(np.int8)
    y = y[:,(((y != 0).sum(axis=0) / y.shape[0]) >= cutoff)]
    return y.astype(np.int8)


def main(args):
    (expr_fn, out_fn) = args[:2]
    expr = np.load(expr_fn)['expr']
    expr1 = nearest_neighbours(expr)
    np.savez_compressed(out_fn, expr=expr1)


if __name__ == '__main__':
    main(sys.argv[1:])
