"""
Download mod l Galois representations from the LMFDB
with image contained in GL_2(F_l).
"""
from sys import argv

if len(argv) != 2:
    print('usage:', argv[0], 'L')
    exit(1)

from lmf import db

query = {'dimension': 2, 'base_ring_characteristic': int(argv[1])}
projection = {'label': True, 'image_type': True,
              'kernel_polynomial': True, 'traces': True}

for x in db.modlgal_reps.search(query, projection):
    print(x)
