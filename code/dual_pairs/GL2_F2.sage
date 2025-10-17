from dual_pairs import DualPair, FiniteFlatAlgebra

R.<x> = QQ[]

pari.read('GL2_F2.gp')
dual_pair_GL2_F2 = pari('dual_pair_GL2_F2')

def dual_pair(h):
    [F, G, Phi] = dual_pair_GL2_F2(h)
    A = FiniteFlatAlgebra(QQ, [R(f) for f in F])
    B = FiniteFlatAlgebra(QQ, [R(g) for g in G])
    return DualPair(A, B, Phi.sage())

def match():
    with open('GL2_F2-reps.txt') as f:
        reps = [eval(x) for x in f.readlines()]

    for r in reps:
        D = dual_pair(R(r['kernel_polynomial']))
        assert D.nice_model().frobenius_traces_lmfdb() == r['traces']
        print(r['label'] + ': ' + str(D.lmfdb_data()))
