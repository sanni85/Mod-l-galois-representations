"""
Compute Frobenius traces of dual pairs
"""
from glob import glob
from sage.arith.misc import primes
from dual_pairs.dual_pair_import import dual_pair_import

def traces(D, B=100):
    # variant of D.frobenius_traces(B)
    P = D.ramified_primes()
    return [-1 if p in P else D.frobenius_matrix(p).trace()
            for p in primes(B)]

def compute_traces(l):
    filenames = glob(str(l^2 - 1) + '.*.gp')
    for f in filenames:
        label = f.replace('.gp', '')
        t = traces(dual_pair_import(f).nice_model())
        print('"' + label + '"\t' + str(t))
