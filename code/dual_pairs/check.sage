"""
Check traces of dual pairs
"""
from dual_pairs import DualPair, FiniteFlatAlgebra

def check(s):
    R = PolynomialRing(QQ, 'x')

    with open(s + '-reps.tsv') as f:
        reps = (x.removesuffix('\n').split('\t') for x in f)
        tr = {label: eval(traces) for label, traces, _, _ in reps}

    with open(s + '-dual-pairs.tsv') as f:
        pairs = (x.removesuffix('\n').split('\t') for x in f)
        for label, data in pairs:
            F, G, d_Phi = eval(data)
            A = FiniteFlatAlgebra(QQ, [R(f) for f in F])
            B = FiniteFlatAlgebra(QQ, [R(g) for g in G])
            Phi = matrix(d_Phi[1]) / d_Phi[0]
            D = DualPair(A, B, Phi)
            if D.nice_model().frobenius_traces_lmfdb() != tr[label]:
                raise ValueError('mismatch for label ' + label)
