# Code to match mod-ell galois representation labels with newform labels
#
# Files mod2dim2reps.txt and mod3dim2reps.txt contain data downloaded
# from the database (table modlgal_reps) using
#
# sage: from lmfdb import db
# sage: modlgal = db.modlgal_reps
# sage: cols = ['label', 'good_primes', 'frobenius_matrices', 'kernel_polynomial', 'projective_kernel_polynomial']
# sage: modlgal.copy_to("mod2dim2reps.txt", columns=cols, include_id=False, query={'base_ring_characteristic':2, 'dimension':2})
# sage: modlgal.copy_to("mod3dim2reps.txt", columns=cols, include_id=False, query={'base_ring_characteristic':3, 'dimension':2})
#
# As of 2024-08-21 this resulted in 503 lines for ell=2 and 6880 lines
# for ell=3, excluding 3 header lines in each output file.  The columns are output in the order
# frobenius_matrices|good_primes|kernel_polynomial|label|projective_kernel_polynomial
#

def read_database_dump(filename):
    data = {}
    n = 0
    with open(filename) as input:
        for L in input:
            n += 1
            if n<4:
                continue
            frob_mats, good_primes, kpol, label, pkpol = L.split("|")
            ell, d, N, i = [ZZ(x) for x in label.split(".")]
            F = GF(ell)
            frob_mats = frob_mats[2:-2].split('},{')
            frob_mats = [[F(a) for a in m.split(",")] for m in frob_mats]
            traces = [m[0]+m[3] for m in frob_mats]
            good_primes = [ZZ(p) for p in good_primes[1:-1].split(",")]
            ap = dict(zip(good_primes, traces))
            kpol = [ZZ(p) for p in kpol.replace("{","").replace("}","").split(",")]
            pkpol = [ZZ(p) for p in pkpol.replace("{","").replace("}","").split(",")]
            data[label] = {'ell':ell, 'N': N, 'ind':i, 'label':label, 'ap': ap, 'kpol':kpol, 'pkpol':pkpol}
    return data

def read_mfmod_data(filename, nap=None):
    data = {}
    label_counter = {}
    with open(filename) as input:
        for L in input:
            label, N, k, c, ind, dim, ell, index, chi_mod_ell, ap = L.split(":")
            # The newform label is N.k.c.ind, but each one can have more than one mod-ell reduction!
            # So we must decorate the label of the reductions to make them unique
            ell = ZZ(ell)
            N = ZZ(N)
            ap = [ZZ(a) for a in ap.split(",")]
            if nap and nap<len(ap):
                ap = ap[:nap]
            ap = dict(zip(primes_first_n(len(ap)), ap))
            i = label_counter.get(label, 0) + 1
            label_counter[label] = i
            xlabel = "-".join([label, str(i)])
            data[xlabel] = {'ell':ell, 'label':label, 'xlabel':xlabel, 'ap': ap}
    # later we'll also read the kpol and pkpol from the files mod*bb*
    return data

def compare(data1, data2):
    """
    Return True if
    (1) data1['ell']==data2['ell']
    (2) for all p in both data1['ap'].keys() and data2['ap'.keys()], data1['ap'][p]==data2['ap'][p]
    """
    if data1['ell'] != data2['ell']:
        return False
    for p, ap in data1['ap'].items():
        if p in data2['ap']:
            if ap != data2['ap'][p]:
                return False
    return True

def match_labels(file1, file2, outputfile, reverse=False):
    """
    file1 contains a database dump for some ell
    file2 contains the modular form reductions modulo ell

    If not reverse: outputfile will contain one line for each label in
    file1. Each line has the label from file1, then ':', then a
    comma-separated list of labels from file 2 which match it.

    If reverse: outputfile will contain one line for each label in
    file2. Each line has the label from file2, then ':', then a
    comma-separated list of labels from file 1 which match it.

    """
    data1 = read_database_dump(file1)
    data2 = read_mfmod_data(file2)
    if reverse:
        data1, data2 = data2, data1
    n = 0
    with open(outputfile, 'w') as out:
        for l1, m1 in data1.items():
            matches = [l2 for l2, m2 in data2.items() if compare(m1,m2)]
            out.write(f"{l1}:{','.join(matches)}\n")
            n += 1
    print(f"{n} lines written to {outputfile}")
