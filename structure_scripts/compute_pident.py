#!/usr/bin/env python
from Bio import SeqIO
import itertools
import numpy as np
import pandas as pd
import sys

in_fn = sys.argv[1]

recs = {rec.id: np.frombuffer(str(rec.seq).encode('ascii'), 'S1')
        for rec in SeqIO.parse(in_fn, 'fasta')}

def gen():
    for (g1, g2) in itertools.combinations(recs, 2):
        s1 = recs[g1]
        s2 = recs[g2]
        f = ((s1 == b'-') & (s2 == b'-'))
        s1 = s1[~f]
        s2 = s2[~f]
        aln_len = len(s1)
        match = (s1 == s2).sum()

        # compute % ident
        yield (g1, g2, 100*(match/aln_len))

df = pd.DataFrame(gen(), columns=['g1', 'g2', 'pident'])
df.to_csv(sys.stdout, sep='\t', index=False)

