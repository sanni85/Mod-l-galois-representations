"""
Compute Frobenius traces of dual pairs
"""
from sage.arith.misc import primes
from dual_pairs.dual_pair_import import dual_pair_import

def traces(D, B=100):
    # variant of D.frobenius_traces(B)
    P = D.ramified_primes()
    return [-1 if p in P else D.frobenius_matrix(p).trace()
            for p in primes(B)]

def compute_traces(f):
    print(traces(dual_pair_import(f).nice_model()))
