"""
Check traces of dual pairs
"""
from dual_pairs import DualPair, FiniteFlatAlgebra

def check(s):
    R = PolynomialRing(QQ, 'x')

    with open(s + '-reps.tsv') as f:
        reps = (x.split('\t') for x in f.readlines())
        tr = {r[0]: eval(r[1]) for r in reps}

    with open(s + '-dual-pairs.tsv') as f:
        for x in f.readlines():
            m = x.split('\t')
            label = m[0]
            F, G, d_Phi = eval(m[1])
            A = FiniteFlatAlgebra(QQ, [R(f) for f in F])
            B = FiniteFlatAlgebra(QQ, [R(g) for g in G])
            Phi = matrix(d_Phi[1]) / d_Phi[0]
            D = DualPair(A, B, Phi)
            if D.nice_model().frobenius_traces_lmfdb() != tr[label]:
                raise ValueError('mismatch for label ' + label)
