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
    # Use to read the mf reduction data from the files 1-1000.txt & similar
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
    return data

ZZx = PolynomialRing(ZZ, 'x')

def parse_poly_string(pol):
    """
    Given a polynomial in ZZ[x], return the list of its coefficients.
    """
    return ZZx(pol).coefficients(sparse=False)

def read_mfmod_bbdata(filename, ell, ignore_reducibles=True):
    """
    Use to read the kpol and pkpol from the files mod[23]bb*.  If
    ignore_reducibles is True, only the irreducible ones are returned.

    For ell=2 the format of each line is one of

    label:reducible
    label:image_name:kernel_poly

    where image_name is 'S3' or 'C3' and 'kernel_poly' is a polynomial in x as a string e.g. 'x^3-x^2-x-1'.

    For ell=3 the format of each line is one of

    label:reducible
    label:image_name:kernel_poly:projective_image_name:projective_kernel_poly

    where image_name is 'GL(2,3)' or 'Nn' or 'Ns', projective_image_name is (repectively) 'S4', D4, V4';
    and 'kernel_poly', 'projective_kernel_poly' are polynomial strings e.g. 'x^3-x^2-x-1'.

    """
    if ell not in [2,3]:
        raise RuntimError(f"BBdata only for ell=2, 3, not {ell}")
    data = {}
    with open(filename) as input:
        for L in input:
            L = L.strip()
            line_data = L.split(":")
            label = line_data[0]
            image = line_data[1]
            if ell==2:
                if len(line_data)>2:
                    kpol = parse_poly_string(line_data[2])
                else:
                    kpol = None
                if not (image == 'reducible' and ignore_reducibles):
                    data[label] = {'ell':ell, 'label':label,
                                   'image':image, 'kpol': kpol}
            if ell==3:
                if len(line_data)>2:
                    kpol = parse_poly_string(line_data[2])
                    proj_image = line_data[3]
                    pkpol = parse_poly_string(line_data[4])
                else:
                    kpol = pkpol = proj_image = None
                if not (image == 'reducible' and ignore_reducibles):
                    data[label] = {'ell':ell, 'label':label,
                                   'image':image, 'kpol': kpol,
                                   'proj_image':proj_image, 'pkpol': pkpol}
    return data

def compare_ap(data1, data2):
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

# dict with keys pairs of polynomials (as coefficient sequences), values True iff splitting fields are isomorphic
splitting_field_dict = {}

def equal_splitting_fields(pol1, pol2):
    global splitting_field_dict
    k = (tuple(pol1), tuple(pol2))
    if k in splitting_field_dict:
        return splitting_field_dict[k]
    K1 = ZZx(pol1).splitting_field('a1')
    K2 = ZZx(pol2).splitting_field('a2')
    res = K1.is_isomorphic(K2)
    splitting_field_dict[k] = res
    return res

def compare_fields(data1, data2):
    """
    Return True if
    (1) data1['ell']==data2['ell']
    (2) data1['kpol']==data2['kpol']
        or they have the same splitting field
    (3) [ell=3 only] data1['pkpol']==data2['pkpol']
        or they have the same splitting field

    """
    ell = data1['ell']
    #print(f"comparing {data1} and {data2} which have compatible ap mod {ell}")

    #return ell == data2['ell'] and data1['kpol'] == data2['kpol'] and (ell==2 or data1['pkpol'] == data2['pkpol'])

    if ell!=data1['ell']:
        return False
    if data1['kpol'] != data2['kpol']:
        print(f"{data1['label']} kpol {data1['kpol']} and {data2['label']} kpol {data2['kpol']} differ!")
        print(f"image = {data2['image']} - checking splitting fields")
        if not equal_splitting_fields(data1['kpol'], data2['kpol']):
            print(f"{data1['label']} and {data2['label']}: database and BB have non-isomorphic splitting fields")
            return False
        else:
            print("OK, splitting fields agree")
    if ell==2:
        return True
    if data1['pkpol'] != data2['pkpol']:
        print(f"{data1['label']} pkpol {data1['pkpol']} and {data1['label']} pkpol {data2['pkpol']} differ!")
        print(f"proj.image = {data2['proj_image']} - checking splitting fields")
        if not equal_splitting_fields(data1['pkpol'], data2['pkpol']):
            print(f"{data1['label']} and {data2['label']}: database and BB have non-isomorphic projective splitting fields")
            return False
        else:
            print("OK, projective splitting fields agree")
    return True

def sort_key_modlgal(label):
    """
    Label is ell.d.N.n with all fields numeric
    """
    return [ZZ(x) for x in label.split(".")]

def sort_key_modform(label):
    """
    Label is N.k.x.y with N,k numeric, x and y alphabetic
    """
    t = label.split(".")
    for i in range(2):
        t[i] = ZZ(t[i])
    return t

def sort_key_modformx(label):
    """
    Label is N.k.x.y-i with N,k,i numeric, x and y alphabetic
    """
    t = label.split(".")
    t = t[:3] + t[3].split("-")
    for i in [0,1,4]:
        t[i] = ZZ(t[i])
    return t

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
    outdata = {}
    for l1, m1 in data1.items():
        matches = [l2 for l2, m2 in data2.items() if compare_ap(m1,m2)]
        matches.sort(key = sort_key_modlgal if reverse else sort_key_modformx)
        outdata[l1] = matches

    keys = list(outdata.keys())
    keys.sort(key = sort_key_modformx if reverse else sort_key_modlgal)
    with open(outputfile, 'w') as out:
        for k in keys:
            out.write(f"{k}:{','.join(outdata[k])}\n")
            n += 1

    print(f"{n} lines written to {outputfile}")

def match_labels_and_fields(ell, file1, file2, file3, outputfile):
    """
    ell is 2 or 3
    file1 contains a database dump for some ell
    file2 contains the modular form reductions modulo ell
    file3 contains the BB data from the modular form reductions modulo ell

    Outputfile will contain one line for each label in
    file1. Each line has the label from file1, then ':', then a
    comma-separated list of labels from file 2 which match it.

    """
    if ell not in [2,3]:
        raise RuntimError(f"BBdata only for ell=2, 3, not {ell}")
    data1 = read_database_dump(file1)
    print(f"Read {len(data1)} items from {file1}")
    data2 = read_mfmod_data(file2)
    print(f"Read {len(data2)} items from {file2}")
    data3 = read_mfmod_bbdata(file3, ell) # ignores reducibles
    print(f"Read {len(data3)} items from {file3}")
    n = 0
    outdata = {}
    for l1, m1 in data1.items():
        matches = [l2 for l2, m2 in data2.items() if compare_ap(m1,m2)]
        field_matches = [l2 for l2 in matches if l2 in data3 and compare_fields(m1,data3[l2])]
        matches.sort(key = sort_key_modformx)
        field_matches.sort(key = sort_key_modformx)
        # if matches:
        #     print(f"{l1} matches {matches} on ap")
        #     print(f"{l1} matches {field_matches} on splitting fields")
        # else:
        #     print(f"{l1} matches nothing on ap")
        if matches != field_matches:
            print(f"{l1} matches {matches} on ap but only {field_matches} on splitting fields")
        outdata[l1] = field_matches

    keys = list(outdata.keys())
    keys.sort(key = sort_key_modlgal)
    with open(outputfile, 'w') as out:
        for k in keys:
            out.write(f"{k}:{','.join(outdata[k])}\n")
            n += 1

    print(f"{n} lines written to {outputfile}")
