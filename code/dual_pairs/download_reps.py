"""
Download mod l Galois representations from the LMFDB.
"""
from sys import argv

argc = len(argv)
if argc < 2 or argc > 3:
    print('usage:', argv[0], 'L [N=2]')
    exit(1)

from lmf import db

query = {'dimension': int(argv[2]) if argc > 2 else 2,
         'base_ring_characteristic': int(argv[1])}
projection = {'label': True, 'image_type': True,
              'kernel_polynomial': True, 'traces': True}

for x in db.modlgal_reps.search(query, projection):
    print(x)
