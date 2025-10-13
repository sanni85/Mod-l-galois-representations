from dual_pairs import DualPair, FiniteFlatAlgebra

R.<x> = QQ[]

# copied from traces.sage
def traces(D, B=100):
    # variant of D.frobenius_traces(B)
    P = D.ramified_primes()
    return [-1 if p in P else D.frobenius_matrix(p).trace()
            for p in primes(B)]

pari.default('realprecision', 120)
pari.read('GSp4_F2.gp')
dual_pair_GSp4_F2 = pari('dual_pair_GSp4_F2')
twin = pari('twin')

def compute_dual_pair(h):
    [F, G, Phi] = dual_pair_GSp4_F2(h)
    A = FiniteFlatAlgebra(QQ, [R(f) for f in F])
    B = FiniteFlatAlgebra(QQ, [R(g) for g in G])
    D = DualPair(A, B, Phi.sage())
    return D, traces(D.nice_model())

def dual_pair(h, tr):
    D, t = compute_dual_pair(h)
    if t == tr:
        return D
    D, t = compute_dual_pair(twin(h))
    if t == tr:
        return D
    raise NotImplementedError

# load dual pairs that have already been computed
import re
r = re.compile('(.*): (.*)')
pairs = {}
try:
    with open('GSp4_F2_matches.txt') as f:
        for x in f.readlines():
            m = r.match(x)
            pairs[m[1]] = eval(m[2])
except FileNotFoundError:
    pass

with open('GSp4_F2_reps.txt') as reps:
    for x in reps.readlines():
        r = eval(x)
        l = r['label']
        if l not in pairs:
            D = dual_pair(R(r['kernel_polynomial']), r['traces'])
            print(l + ': ' + str(D.lmfdb_data()))
