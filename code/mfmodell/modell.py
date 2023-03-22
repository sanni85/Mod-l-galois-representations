# To run this: have an lmfdb directory containing the code from
# https://github.com/LMFDB/lmfdb with appropriate configuration
# settings, put this file in there, start Sage and type
#
# "%runfile  modell_new.py"
#
# I only run this from legendre.mit.edu or speed of access to the
# database.
#

from sage.all import ZZ, QQ, GF, PolynomialRing, NumberField, primes, primes_first_n, Matrix, vector, FiniteDimensionalAlgebra
from sage.databases.cremona import class_to_int
from mf_compare import read_dtp
from lmfdb import db
from lmfdb.classical_modular_forms.web_newform import WebNewform
import logging
from lmfdb.classical_modular_forms import cmf_logger
cmf_logger.setLevel(logging.NOTSET)

# Function to check that the order m of some Dirichlet character chi
# satisfies the necessary condition for a newform with character chi
# to have a mod-ell reduction.

def char_order_valid(m, ell):
    """Return True iff the positive integer m has the form (ell^k)*m1
    where m1 is 1 mod ell
    """
    m1 = ZZ(m).prime_to_m_part(ell)
    return m1.divides(ell-1)

# Function to find primes above ell in the Hecke field.  No longer
# used.

def Hecke_field_primes(fpol, ell, degree=None, verbose=False):
    """Given a polynomial fpol defining a number field K and a prime ell,
    return the list of primes of K above ell.  Set degree to a
    suitable positive integer to only return primes of that degree,
    e.g. degree=1 will only give primes whose residue field is
    GF(ell).
    """
    if verbose:
        print("Checking field poly {} at ell={}".format(fpol,ell))
    Qx = PolynomialRing(QQ,'x')
    K = NumberField(Qx(fpol), 'a_')
    ellK = K.primes_above(ell, degree=degree)
    if verbose:
        print("{} primes of norm {}".format(len(ellK),ell))
    return ellK

def make_an_decoder(nf):
    """Given a newform nf, return a function which converts vectors of
    coefficients (in QQ) to elements of the Hecke field of nf.

    Not currently used.
    """
    K = nf.hecke_field
    if K==QQ:
        return lambda elt: QQ(elt)
    assert nf.hecke_ring_cyclotomic_generator==0

    if nf.hecke_ring_power_basis:
        a = K.gen()
        betas = [a**i for i in range(K.degree())]
    else: # generic
        betas = [K(c)/d for c,d in zip(nf.hecke_ring_numerators, nf.hecke_ring_denominators)]
    return lambda elt: sum(c*beta for c, beta in zip(elt, betas))
#
# Notation:
#
# K is a number field (the Hecke field of a newform)
# betas = ZZ-module generators of the Hecke order
# ell = the prime

# Function to construct all homomorphisms from the Hecke order to
# GF(ell), without using the number field at all.  Note that asking
# for the primes of K above ell involves computing the maximal order,
# which can be expensive, even though theoretically we only need an
# ell-maximal order.

def reduction_maps(betas, ell, verbose=False):
    """Input: 

    betas: a list of d ZZ-module generators of an order O in a number
    field K of degree d,

    ell: a rational prime ell,

    Output: a (possibly empty) list of maps ZZ^d = GF(ell) encoding all
    possible ring homomorphisms from O to GF(ell).  Each map is
    encoded as a list of d elements of GF(ell) giving the images of
    the order's basis (the betas).

    Method: convert each beta into a dxd integer matrix giving the
    multiplication-by-beta map from O to itself.  Reduce these mod
    ell, and use them to define a FiniteDimensionalAlgebra A over
    GF(ell).  Find the maximal ideals of A whose quotient field is
    GF(ell) (and not some extension of it).  So the ideal has
    dimension d-1 and the required map is the inner product with a
    basis vector for its complement.
    """
    K = betas[0].parent()
    # print("K = {}".format(K))
    # print("beta== 1? {}".format([b==1 for b in betas]))
    # print("beta==-1? {}".format([b==-1 for b in betas]))
    # assert betas[0]==1
    Fl = GF(ell)
    if verbose:
        print("creating Hecke order mod {} from betas".format(ell))
    coords = K.gen().coordinates_in_terms_of_powers()
    U = Matrix([coords(b) for b in betas]).inverse()
    structure = [(Matrix([coords(bi*bj) for bj in betas])*U).change_ring(Fl)
                 for bi in betas]
    one = [ZZ(c) for c in vector(coords(K(1))) * U]

    A = FiniteDimensionalAlgebra(Fl, structure)

    MM = A.maximal_ideals()
    d = len(betas)
    vv = [list(M.basis_matrix().right_kernel().basis()[0])
          for M in MM  if M.basis_matrix().rank()==d-1]

    # Taking the dot product with each of these basis vectors gives a
    # GF(l)-linear map from the Hecke order to GF(l), but we need it
    # to be a ring homomorphism, so we must scale it so 1 maps to 1.
    # In practice Sage will always scale the basis vector so that the
    # first nonzero entry (which here must be the first entry) is 1,
    # but we should not rely on this.  We also do not rely on any of
    # the basis elements (betas) beinb 1.

    for v in vv:
        image_of_one = sum([r*s for r,s in zip(one,v)])
        if image_of_one==0:
            raise ValueError("reduction map defined by {} maps 1 to 0".format(v))
        if image_of_one!=1:
            if verbose:
                print("Rescaling reduction vector v={}".format(v))
            inv = 1/image_of_one
            v = [vi*inv for vi in v]
            if verbose:
                print("Rescaled reduction vector v={}".format(v))
    if verbose:
        print("{} reductions found from Hecke order".format(len(vv)))
        if vv:
            print("Reduction vector(s):")
            for v in vv: print(v)
    return vv

# Function to construct all the reuductions of one newform at a prime
# ell:

def make_reductions(nf, ell, verbose=False):
    """Input: A newform nf and a prime ell

    Output: a list of all the associated reducton maps

    Method: depends on the Hecke field.  Trivial if it is QQ,
    otherwise construct the Hecke order's ZZ-basis (betas) depending
    on the type of basis stored in the newform.

    """
    K = nf.hecke_field
    if verbose:
        print("Making reduction mod ell={}".format(ell))
    if K==QQ:
        nf.betas = [[1]] # for consistency, not actually used
        Fl = GF(ell)
        return [Fl(1)]
    
    try:
        betas = nf.betas # we'll have them already if the data came from a file
    except AttributeError:
        if not nf.hecke_ring_numerators:#nf.hecke_ring_power_basis:
            a = K.gen()
            betas = [a**i for i in range(K.degree())]
        else: # generic
            betas = [K(c)/d for c,d in zip(nf.hecke_ring_numerators, nf.hecke_ring_denominators)]
        nf.betas = betas
    return reduction_maps(betas, ell, verbose)

def apply_red(coeffs, red):
    """Apply a reduction map to one Hecke ring element

    Input:

    - coeffs: list of d integers (d>0) giving the ZZ-coefficients of
      an element of the Hecke ring.

    - red: list of d elements of GF(l) determining a ring homomorphism
      from the Hecke ring to GF(l).

    Output: an element of GF(l).
    """
    # When each red was implemented as a map we would just need to do
    # return red(coeffs)
    
    Fl = red[0].parent()    
    return Fl(sum([ci*ri for ci,ri in zip(coeffs,red)]))

def reduce_an_mod_ell(nf,red):
    """Input: a newform and one reduction map to GF(ell) for some ell

    Output: a list of the reductions (a_n mod ell) of all the
    q-expansion coefficients a_n of the newform.

    """
    return [apply_red(an,red) for an in nf.an]
    
def reduce_ap_mod_ell(nf, red, verbose=False):
    """Input: a newform and one reduction map to GF(ell) for some ell

    Output: a list of pairs (p,a_p mod ell) for the prime index
    q-expansion coefficients a_p of the newform.

    """
    if verbose:
        print("Reducing {} mod ell via {}".format(nf.label, red))
        for i,p in enumerate(primes_first_n(10)):
            ap = nf.ap[i]
            print("p={}, ap={}, ap mod ell = {}".format(p, ap, apply_red(ap, red)))
    return dict([(p,apply_red(nf.ap[i], red)) for i,p in enumerate(primes_first_n(len(nf.ap)))])

def reduce_chi_mod_ell(nf, red):
    """Input: a newform and one reduction map to GF(ell) for some ell

    -  If nf has nontrivial character chi then
    nf.hecke_ring_character_values is a list of pairs [a,coeffs] where
    the a's generate (Z/NZ)^* and the coeffs are the coefficients of
    chi(a) with respect to the Hecke ring basis.

    - If chi is trivial then nf.hecke_ring_character_values==None.

    Output: a list of pairs [a, chi(a) mod ell] for generators a as
    above (empty list for trivial chi).

    """
    if nf.char_order==1:
        return []

    return [[a,apply_red(coeffs, red)] for a, coeffs in nf.hecke_ring_character_values]

########################################################################
#
# Main function get_forms(N,k,ell) finds all mod-ell reductions of
# forms of level N and weight k for which the Hecke field information
# is in the database.
#
########################################################################

def get_forms(N, k, ell, max_dim=50, verbose=False):
    """
    Input: N (level), k (weight), ell (prime); and a verbosity flag

    Output: 2 lists of newforms of level N, weight k, whose character
    chi has order compatible with ell.

    (1) forms for which we have the Hecke field and which have at
    least one mod-ell reduction.  The reductions are stored as an
    attribute (.reductions) for each of these newforms.

    (2) forms for which we do not have the Hecke field in the
    database, for which we currently do nothing.

    """
    if N%ell==0:
        return [], []
    nfs  = db.mf_newforms
    heckes = db.mf_hecke_nf # for more ap
    forms  = [WebNewform.by_label(f_label) for f_label in
              nfs.search({'level':N, 'weight':k}, projection='label')]
    if verbose:
        print("forms with (N,k)=({},{}): {}".format(N,k,[f.label for f in forms]))

    # omit forms whose character order is invalid:
    
    forms = [f for f in forms if char_order_valid(f.char_order, ell)]
    if verbose:
        print("After char order check, forms with (N,k,ell)=({},{},{}): {}".format(N,k,ell,[(f.label,f.dim) for f in forms]))

    # for forms with no Hecke field, try to read from a data file

    forms_with_no_field = [f for f in forms if f.field_poly is None]

    if forms_with_no_field:
        if verbose:
            print("Before reading data files, {} forms have no Hecke field data: {}".format(len(forms_with_no_field), [f.label for f in forms_with_no_field]))

        for f in forms:
            if f.field_poly is None and f.dim <= max_dim:
                if verbose:
                    print("looking for a data file for form {}".format(f.label))
                f = get_form_data_from_file(f, max_dim=max_dim, verbose=verbose)

        forms_with_no_field = [(f.label,f.dim,ell) for f in forms if f.field_poly is None]
        if verbose:
            print("After reading data files, {} forms have no Hecke field data".format(len(forms_with_no_field)))

    # Now exclude any forms which still have no Hecke data:

    forms = [f for f in forms if f.field_poly]
    
    # find all mod-ell reduction maps from the Hecke order

    Qx = PolynomialRing(QQ,'x')
    for f in forms:
        f.hecke_field = NumberField(Qx(f.field_poly), 'a_')
        if verbose:
            print("making reductions for {} mod {}".format(f.label, ell))
            print("Hecke field is {}".format(f.hecke_field))
        f.reductions = make_reductions(f, ell, verbose)
        if verbose:
            print("finished making reductions for {} mod {}".format(f.label, ell))

        try:
            assert f.an
        except AttributeError:
            # get extra ap from the second table:
            anap = heckes.lucky({'label':f.label}, projection=['ap','an'])
            f.an = anap['an']
            f.ap = anap['ap']

            if verbose:
                nap = len(f.ap)
                print("Found {} ap and {} an in the second table".format(nap,len(f.an)))
        
    # Exclude forms with no mod-ell reductions:
        
    forms = [f for f in forms if f.reductions]
    if forms and verbose:
        print("{}.{} forms with mod-{} reductions:".format(N,k,ell))
        for f in forms:
            print("{} has {} reductions mod {}".format(f.label,len(f.reductions), ell))

    return forms, forms_with_no_field

########################################################################
#
# Function to get a single form from a data file output by the Magma
# OneSpace function in bigspaces.m, which just calls NewspaceData and
# outputs to a file.
#
# The filename is expected to be a space label N.k.c or a newform
# label N.k.c.i.  In either case it contains data for all newforms in
# the space N.k.c and we use i (converted to an integer) to know which
# newform we want from that space.
#
########################################################################

DATA_DIR="/scratch/home/jcremona/Mod-l-galois-representations/data/mfmodell"
NF_DIR="/scratch/home/jcremona/Mod-l-galois-representations/data/mfmodell/0"

def get_form_data_from_file(nf, data_dir=NF_DIR, max_dim=50, verbose=False):
    """Fill in data for an incomplete WebNewform object by reading from a
    data file.  If a suitable file does not exist, the original
    WebNewform object is returned unchanged, with a message output.

    The data files only contain the extra data for dimensions less
    than some bound.  If the newform's dimension is greater then we do
    nothing.

    """
    if not nf.field_poly is None:
        return nf
    fname = label = nf.label
    datafile = "/".join([data_dir, fname])
    try:
        data = read_dtp(datafile, verbose=False)
    except FileNotFoundError:
        fname = fname[:fname.rindex(".")]
        datafile = "/".join([data_dir, fname])
        try:
            data = read_dtp(datafile, verbose=False)
        except FileNotFoundError:
            print("No file {} or {} found: no data available".format(datafile,label))
            return nf

    if verbose:
        print("Successfully read data for {} from file {}".format(label, fname))
        
    # Now we have data.  It is a dict with a single key (N,k,o) where
    # o is the character orbit number and value so we just extract the
    # single value for this key, which is another dict:

    data = list(data.values())[0]

    # This dict has keys ['dims', 'traces', 'ALs', 'polys',
    # 'eigdata'], and each value is a list, one per newform, so we
    # need to extract just one item from the relevant lists, according
    # to the newform we want.  This is determined by the 4th (last)
    # part of the label we have, which we need to convert to an
    # integer (from 0).

    # NB eigdata is a list of dicts (see later comment for details)
    # but *only* for components of dimension>1
    
    nf_index = class_to_int(label.split(".")[3])
    dims = data['dims']
    nf_eigdata_index = nf_index - dims.count(1)
    if verbose: # debug
        print("Newform label = {}".format(label))
        print("nf_index = {}".format(nf_index))
        print("len(polys) = {}".format(len(data['polys'])))
        print("len(dims) = {}".format(len(dims)))
        print("#dims>1 = {}".format(len(dims)-dims.count(1)))
        print("len(eigdata) = {}".format(len(data['eigdata'])))
        print("nf_eigdata_index = {}".format(nf_eigdata_index))
    dim = data['dims'][nf_index]
    assert dim == nf.dim
    nf.field_poly = data['polys'][nf_index]

    eigdata = data['eigdata'][nf_eigdata_index]

    # eigdata is a dict with keys  ['poly', 'basis', 'n', 'm', 'ans', 'char']
    chi_gens, chi_vals = eigdata['char']
    nf.hecke_ring_character_values = list(zip(chi_gens, chi_vals))
    # we will not need to access the char_order since this space has
    # already been selected to have an appropriate character.

    if verbose:
        print("hecke_ring_character_values = {}".format(nf.hecke_ring_character_values))
    
    Qx = PolynomialRing(QQ,'x')
    nf.hecke_field = K = NumberField(Qx(nf.field_poly), 'a_')
    nf.betas = [K(b) for b in eigdata['basis']]
    nf.an = eigdata['ans']
    nf.ap = [nf.an[p-1] for p in primes(len(nf.an))]

    if verbose: # debug
        print("We have {} a_n and {} a_p for newforms {}".format(len(nf.an), len(nf.ap), label))
    return nf

########################################################################
#
# Some utility functions
#
########################################################################

def nf_mod_ell(nf, red):
    """Input: a newform and one mod-ell reduction.

    Output: a newform-mod-ell dict (used for output) with entries for
    the label, level, weight, dimension, reduced ap's and reduced
    character values.

    """
    return {'label': nf.label,
            'N': nf.level,
            'k': nf.weight,
            'd': nf.dim,
            'ap': reduce_ap_mod_ell(nf, red),
            'chi_mod_ell': reduce_chi_mod_ell(nf, red),
            }

def compare_record(nf1,nf2):
    """
    Return True iff the two newform-mod-ell dicts have the same ell and the sam vector of ap mod ell.
    """
    return (nf1['ell'] == nf2['ell']) and (nf1['ap'] == nf2['ap'])

def nf_to_string(nfr, nap=None):
    """Return the formatted output string for a newform mod ell dict.  If
    nap is a positive integer (and not None, the default) then only a
    maximum of nap ap's are output.

    The string uses ":" as separators between the fields

    - N   |
    - k   | (the four components of the newform label:
    - c   |   level,weight,char.orbit, nf id)
    - id  |
    - dim         (dimension of the Hecke field)
    - ell         (the prime ell)
    - index       (index of this mod-ell newform among all with the same ap mod ell)
    - chi_mod_ell (character value list mod ell)
    - ap          (list of ap mod ell)

    """
    label = nfr['label']
    N,k,c,id = label.split(".")
    ell = str(nfr['ell'])

    # The ap mod ell are in a single comma-separated list:
    aplist = [nfr['ap'][p] for p in sorted(nfr['ap'].keys())]
    if nap and nap<len(aplist):
        aplist = aplist[:nap]
    ap = str(aplist).replace(" ","")[1:-1]

    index = str(nfr['index'])
    dim = str(nfr['d'])
    chi_mod_ell = str(nfr['chi_mod_ell']).replace(" ","")
    return ":".join([label, N,k,c,id, dim, ell, index, chi_mod_ell, ap])

def data_output(nf_list, filename, mode='w', nap=None):
    """Output a list of full nf-mod-ell records to the given file.  If
    nap is a positive integer (and not None, the default) then only a
    maximum of nap ap's are output.  NB: contents of file will be
    over-written unless mode='a'.

    """
    o = open(filename, mode=mode)
    for nf in nf_list:
        o.write(nf_to_string(nf, nap) + "\n")
    o.close()
    nuniq = len([nf['label'] for nf in nf_list if nf['index']==1])
    print("{} forms (with {} unique reductions) output to {}".format(len(nf_list),nuniq,filename))

def extra_output(nf_list, filename, mode='w'):
    """Output a list of basic newform info to the given file.

    The only output label:dimension:ell for each, and is used for
    newforms which may have mod-ell reduction but we have not computed
    them, for example if we do not have thir Hecke field.

    """
    o = open(filename, mode=mode)
    for nf in nf_list:
        # each is a tuple (label, dim, ell)
        lab = nf[0]
        space = lab[:lab.rfind(".")]
        id = lab[1+lab.rfind("."):]
        o.write("{}:{}:{}:{}\n".format(space, id, nf[1], nf[2]))
    o.close()
    print("{} unprocessed forms output to {}".format(len(nf_list),filename))

########################################################################
#
# Function for a systematic run over a range of levels and primes
#
# After running, use data_output() and extra_output() to output to suitable files
#
# e.g.
# sage: res = run([1..100],[2,3,5], verbose=True)
# sage: for ell in [2,3,5]:
#           data_output(res[0][ell], "mod_{}_100.txt".format(ell))
#           extra_output(res[1][ell], "mod_{}_100_missing.txt".format(ell))
#
########################################################################
    
def run(levels, ells=[2,3,5], max_dim=50, verbose=True):
    """
    Input:

    - levels: a list of levels
    - ells: a list pr primes (default [2,3,5])
    - verbose: verbosity flag (default True)

    Output:

    tuple of two dicts nf and nfx, each with the ells as keys.

    -nf[ell]  is a list of full nf_mod_ell records
    -nfx[ell] is a list of triple (label,dimension,ell) of unprocessed newforms

    """
    nf = dict([(ell,[]) for ell in ells])
    nnf = dict([(ell,0) for ell in ells])
    nfx = dict([(ell,[]) for ell in ells])

    for ell in ells:
        for N in levels:
            for k in range(2,1+max(ell+1,4)):
                if True:#verbose:
                    print("(ell,k,N) = ({},{},{})".format(ell,k,N))
                ff, ff_no_field = get_forms(N,k,ell, max_dim=max_dim, verbose=verbose)
                if ff_no_field:
                    if verbose:
                        print("No Hecke field exists for {}".format(ff_no_field))
                    nfx[ell] += ff_no_field
                if ff and verbose:
                    print("(ell,k,N) = ({},{},{})".format(ell,k,N))
                for f in ff:
                    for red in f.reductions:
                        nfr = nf_mod_ell(f, red)
                        nfr['ell'] = ell
                        # See how many (if any) times we have seen this list of ap mod ell before
                        n_previous = sum([compare_record(nfr,nfi) for nfi in nf[ell]])
                        if n_previous==0:
                            nnf[ell] += 1 # count distinct mode ell forms
                        nfr['index'] = n_previous+1
                        nf[ell].append(nfr) # collect all, including repeats
                        if verbose:
                            if n_previous==0:
                                print("New mod {} form from {} of dimension {}".format(ell,nfr['label'],nfr['d']))
                            else:
                                print("Repeat (#{}) mod {} form from {} of dimension {}".format(n_previous+1,ell,nfr['label'],nfr['d']))

        print("{} mod {} newforms found, with {} distinct".format(len(nf[ell]), ell, nnf[ell]))
        if nfx[ell]:
            if verbose:
                for f in nfx[ell]:
                    print("Unable to reduce {} mod {} (dimension={})".format(f[0],f[2],f[1]))
            else:
                print("Unable to reduce {} newforms mod {}".format(len(nfx[ell]), ell))
    return nf, nfx

def compare_new_old(new_f, old_f):
    # compare levels
    if old_f['N'] != new_f['N']:
        return False
    for p in new_f['ap_mod_ell'].keys():
        if new_f['ap_mod_ell'][p] != old_f['ap'][p]:
            return False
    return True

"""
To extract the dimensions from an "extra" output file:
awk -F ":" '{print $3;}'  mod_5_100_missing.txt | sort -n | uniq

To count the distinct mod-ell reductions in a normal output file:
awk -F ":" '$8==1{print $8;}' mod_5_100.txt | wc -l

To see the largest multiplicities:
awk -F ":" '$8>1{print $8;}' mod_2_100.txt | sort -n | uniq | tail

"""
