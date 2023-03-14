from __future__ import print_function

import sys

from Bio import SeqIO
from ipython_genutils.py3compat import xrange

rec = SeqIO.read(sys.argv[1], 'genbank')
feats = [feat for feat in rec.features if feat.type == "CDS"]
start, end = sys.argv[2].split(':')

desired = set(xrange(int(start), int(end), 1))

for f in feats:
    span = set(xrange(f.location._start.position, f.location._end.position))
    if span & desired:
        print(f)
