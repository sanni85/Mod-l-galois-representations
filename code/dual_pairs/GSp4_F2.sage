from dual_pairs import DualPair, FiniteFlatAlgebra

R.<x> = QQ[]

pari.default('realprecision', 120)
pari.read('GSp4_F2.gp')
dual_pair_GSp4_F2 = pari('dual_pair_GSp4_F2')
twin = pari('twin')

def compute_dual_pair(h):
    [F, G, Phi] = dual_pair_GSp4_F2(h)
    A = FiniteFlatAlgebra(QQ, [R(f) for f in F])
    B = FiniteFlatAlgebra(QQ, [R(g) for g in G])
    D = DualPair(A, B, Phi.sage())
    return D, D.nice_model().frobenius_traces_lmfdb()

def dual_pair(h, tr):
    D, t = compute_dual_pair(h)
    if t == tr:
        return D
    D, t = compute_dual_pair(twin(h))
    if t == tr:
        return D
    raise NotImplementedError

def match():
    with open('GSp4_F2-reps.txt') as f:
        reps = [eval(x) for x in f.readlines()]

    for r in reps:
        D = dual_pair(R(r['kernel_polynomial']), r['traces'])
        print(r['label'] + ': ' + str(D.lmfdb_data()))
