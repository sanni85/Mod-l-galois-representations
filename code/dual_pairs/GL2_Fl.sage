from dual_pairs.dual_pair_from_sub_GL2_field import dual_pair_from_sub_GL2_field

R.<x> = QQ[]

@cached_function
def candidates(p, F):
    # catch sporadic failure
    for i in range(8):
        try:
            return dual_pair_from_sub_GL2_field(p, F)
        except ValueError:
            continue

def dual_pair_string(F, t, p, image_type, traces_big):
    if image_type == 'big':
        with open(traces_big[tuple(t)] + '.dualpair') as f:
            return f.read().removesuffix('\n')
    Ds = candidates(p, F)
    if len(Ds) == 1:
        return str(Ds[0].lmfdb_data())
    # TODO: pick the correct D using fewer traces
    Ds = [D for D in Ds if D.nice_model().frobenius_traces_lmfdb() == t]
    if len(Ds) != 1:
        raise ValueError('no uniquely defined dual pair')
    return str(Ds[0].lmfdb_data())

def match(l):
    F = GF(l)

    # traces for representations with big image
    with open('GL2_F' + str(l) + '-traces.tsv') as f:
        tr = [s.removesuffix('\n').split('\t') for s in f]
        traces_big = {tuple(eval(y)): x for x, y in tr}
        if len(traces_big) != len(tr):
            raise ValueError('non-unique sequence of traces')

    with open('GL2_F' + str(l) + '-reps.tsv') as f:
        reps = (x.removesuffix('\n').split('\t') for x in f)
        for label, t, p, image_type in reps:
            print(label + '\t' + dual_pair_string(F, eval(t), R(eval(p)),
                                                  image_type, traces_big))
