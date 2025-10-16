"""
Download one-dimensional mod l Galois representations from the LMFDB.
"""
from lmf import db

query = {'dimension': 1}
projection = {'label': True, 'traces': True}

for x in db.modlgal_reps.search(query, projection):
    print(x)
