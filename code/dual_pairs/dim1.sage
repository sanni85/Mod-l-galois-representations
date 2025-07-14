import re

from dual_pairs import DualPair, FiniteFlatAlgebra

abgalrep = pari.read('abgalrep.gp')
r = re.compile('([0-9]+)\\.1\\.[0-9]+\\.([0-9]+)-([0-9]+)')
R.<x> = QQ[]

def data_from_char_label(label):
    m = r.match(label)
    l = int(m[1])
    n = int(m[2])
    c = int(m[3])
    return abgalrep(n, c, l)

def dual_pair_from_data(data):
    F, G, d_Phi = data
    Phi = d_Phi[1] / d_Phi[0]
    A = FiniteFlatAlgebra(QQ, [R(f) for f in F])
    B = FiniteFlatAlgebra(QQ, [R(g) for g in G])
    return DualPair(A, B, Phi.sage())

# copied from traces.sage
def traces(D, B=100):
    # variant of D.frobenius_traces(B)
    P = D.ramified_primes()
    return [-1 if p in P else D.frobenius_matrix(p).trace()
            for p in primes(B)]

def match():
    with open('dim-1-reps.txt') as f:
        reps = [eval(x) for x in f.readlines()]

    for r in reps:
        data = data_from_char_label(r['label'])
        D = dual_pair_from_data(data)
        assert traces(D.nice_model()) == r['traces']
        print(r['label'] + ': ' + str(data))
