from dual_pairs.dual_pair_from_sub_GL2_field import dual_pair_from_sub_GL2_field

R.<x> = QQ[]

# copied from traces.sage
def traces(D, B=100):
    # variant of D.frobenius_traces(B)
    P = D.ramified_primes()
    return [-1 if p in P else D.frobenius_matrix(p).trace()
            for p in primes(B)]

@cached_function
def candidates(p, F):
    return dual_pair_from_sub_GL2_field(p, F)

def dual_pair_string(F, r, traces_big):
    t = r['traces']
    if r['image_type'] == 'big':
        # TODO: replace by grep?
        X = [x for x in traces_big if x[1] == t]
        if len(X) != 1:
            raise ValueError('no uniquely defined dual pair')
        with open(X[0][0] + '.gp') as f:
            s = f.read().removesuffix('\n')
    else:
        p = R(r['kernel_polynomial'])
        Ds = [D for D in candidates(p, F)
              if traces(D.nice_model()) == t]
        if len(Ds) != 1:
            raise ValueError('no uniquely defined dual pair')
        s = str(Ds[0].lmfdb_data())
    return s

def match(l):
    F = GF(l)

    # traces for representations with big image
    with open('mod-' + str(l) + '-traces.txt') as t:
        traces_big = [eval(x) for x in t.readlines()]

    with open('mod-' + str(l) + '-reps.txt') as f:
        reps = [eval(x) for x in f.readlines()]

    for r in reps:
        print(r['label'] + ': ' + dual_pair_string(F, r, traces_big))
