"""
Download number fields from the LMFDB with given Galois group.
"""
from sys import argv

if len(argv) != 2:
    print('usage:', argv[0], 'LABEL')
    exit(1)

from lmf import db

query = {'galois_label': argv[1]}
projection = {'label': True, 'coeffs': True}

for x in db.nf_fields.search(query, projection):
    print('"' + x['label'] + '"\t' + str(x['coeffs']))
