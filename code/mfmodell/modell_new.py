# To run this: have an lmfdb directory containing the code from
# https://github.com/LMFDB/lmfdb with appropriate configuration
# settings, put this file in there, start Sage and type
#
# "%runfile  modell_new.py"
#
# I only run this from legendre.mit.edu or speed of access to the
# database.
#

from sage.all import ZZ, QQ, GF, PolynomialRing, NumberField, primes, Matrix, FiniteDimensionalAlgebra
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

def reduction_maps(betas, ell):
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
    #print("K = {}".format(K))
    #print("betas = {}".format(betas))
    assert betas[0]==1
    Fl = GF(ell)
    print("creating Hecke order mod {} from betas".format(ell))
    coords = K.gen().coordinates_in_terms_of_powers()
    U = Matrix([coords(b) for b in betas]).inverse()
    structure = [(Matrix([coords(bi*bj) for bj in betas])*U).change_ring(Fl)
                 for bi in betas]
    assert structure[0]==1
    #print("structure = {}".format(structure))
    A = FiniteDimensionalAlgebra(Fl, structure)
    #print("A = {}".format(A))
    MM = A.maximal_ideals()
    d = len(betas)
    vv = [list(M.basis_matrix().right_kernel().basis()[0])
          for M in MM  if M.basis_matrix().rank()==d-1]
    print("{} reductions found from Hecke order".format(len(vv)))
    return [lambda w: Fl(sum([vi*Fl(wi) for vi,wi in zip(v,w)])) for v in vv]

# Function to construct all the reuductions of one newform at a prime
# ell:

def make_reductions(nf, ell):
    """Input: A newform nf and a prime ell

    Output: a list of all the associated reducton maps

    Method: depends on the Hecke field.  Trivial if it is QQ,
    otherwise construct the Hecke order's ZZ-basis (betas) depending
    on the type of basis stored in the newform.

    """
    K = nf.hecke_field
    if K==QQ:
        Fl = GF(ell)
        return [lambda elt: Fl(QQ(elt))]
    
    assert nf.hecke_ring_cyclotomic_generator==0

    if nf.hecke_ring_power_basis:
        a = K.gen()
        betas = [a**i for i in range(K.degree())]
    else: # generic
        betas = [K(c)/d for c,d in zip(nf.hecke_ring_numerators, nf.hecke_ring_denominators)]
    return reduction_maps(betas, ell)
   
def reduce_an_mod_ell(nf,red):
    """Input: a newform and one reduction map to GF(ell) for some ell

    Output: a list of the reductions (a_n mod ell) of all the
    q-expansion coefficients a_n of the newform.

    """
    return [red(an) for an in nf.qexp]
    
def reduce_ap_mod_ell(nf, red):
    """Input: a newform and one reduction map to GF(ell) for some ell

    Output: a list of pairs (p,a_p mod ell) for the prime index
    q-expansion coefficients a_p of the newform.

    """
    return dict([(p,red(nf.qexp[p])) for p in primes(nf.qexp_prec)])

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
    return [[a,red(coeffs)] for a, coeffs in nf.hecke_ring_character_values]

########################################################################
#
# Main function get_forms(N,k,ell) finds all mod-ell reductions of
# forms of level N and weight k for which the Hecke field information
# is in the database.
#
########################################################################

def get_forms(N, k, ell, verbose=False):
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
    forms  = [WebNewform.by_label(f_label) for f_label in
              nfs.search({'level':N, 'weight':k}, projection='label')]
    if verbose:
        print("forms with (N,k)=({},{}): {}".format(N,k,[f.label for f in forms]))

    # omit forms whose character order is invalid:
    
    forms = [f for f in forms if char_order_valid(f.char_order, ell)]
    if verbose:
        print("After char order check, forms with (N,k,ell)=({},{},{}): {}".format(N,k,ell,[(f.label,f.dim) for f in forms]))

    # identify forms with no Hecke field:
    
    forms_with_no_field = [(f.label,f.dim,ell) for f in forms if f.field_poly==None]

    # Now exclude these:

    forms = [f for f in forms if f.field_poly]
    
    # find all mod-ell reduction maps from the Hecke order

    Qx = PolynomialRing(QQ,'x')
    for f in forms:
        f.hecke_field = NumberField(Qx(f.field_poly), 'a_')
        if verbose:
            print("making reductions for {} mod {}".format(f.label, ell))
        f.reductions = make_reductions(f, ell)
        if verbose:
            print("finished making reductions for {} mod {}".format(f.label, ell))

    # Exclude forms with no mod-ell reductions:
        
    forms = [f for f in forms if f.reductions]
    if forms and verbose:
        print("{}.{} forms with mod-{} reductions:".format(N,k,ell))
        for f in forms:
            print("{} has {} reductions mod {}".format(f.label,len(f.reductions), ell))

    return forms, forms_with_no_field

########################################################################
#
# Some utility functions
#
########################################################################

def nf_mod_ell(nf, red):
    """
    Input: a newforms and a mod-ell reduction

    Output: a newform-mod-elll dictionary (used for output) with
    entries for the label, level, weight, dimension, reduced ap's and
    reduced character values.

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

def nf_to_string(nfr):
    """
    Return the formatted output string for a newform mod ell dict.

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
    # old long format
    # ap = ["[{},{}]".format(p,nfr['ap'][p]) for p in sorted(nfr['ap'].keys())]

    # The ap mod ell are in a single comma-separated list:
    ap = str([nfr['ap'][p] for p in sorted(nfr['ap'].keys())]).replace(" ","")[1:-1]

    index = str(nfr['index'])
    dim = str(nfr['d'])
    chi_mod_ell = str(nfr['chi_mod_ell']).replace(" ","")
    return ":".join([label, N,k,c,id, dim, ell, index, chi_mod_ell, ap])

def data_output(nf_list, filename, mode='w'):
    """
    Output a list of full nf-mod-ell records to the given file.

    NB: contents of file will be over-written unless mode='a'.
    """
    o = open(filename, mode=mode)
    for nf in nf_list:
        o.write(nf_to_string(nf) + "\n")
    o.close()
    print("{} forms output to {}".format(len(nf_list),filename))

def extra_output(nf_list, filename, mode='w'):
    """Output a list of basic newform info to the given file.

    The only output label:dimension:ell for each, and is used for
    newforms which may have mod-ell reduction but we have not computed
    them, for example if we do not have thir Hecke field.

    """
    o = open(filename, mode=mode)
    for nf in nf_list:
        # each is a tuple (label, dim, ell)
        o.write("{}:{}:{}\n".format(nf[0],nf[1],nf[2]))
    o.close()
    print("{} unprocessed forms output to {}".format(len(nf_list),filename))

########################################################################
#
# Function for a systematic run over a range of levels and primes
#
########################################################################
    
def run(levels, ells=[2,3,5], verbose=True):
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
        for k in range(2,1+max(ell+1,4)):
            for N in levels:
                if verbose:
                    print("(ell,k,N) = ({},{},{})".format(ell,k,N))
                ff, ff_no_field = get_forms(N,k,ell)
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
            for f in nfx[ell]:
                print("Unable to reduce {} mod {} (dimension={})".format(f[0],f[2],f[1]))
    return nf, nfx

# After running the above use data_output() to output nf[ell] for each ell to a suitable file, and extra_output() to output nfx[ell]

def compare_new_old(new_f, old_f):
    # compare levels
    if old_f['N'] != new_f['N']:
        return False
    for p in new_f['ap_mod_ell'].keys():
        if new_f['ap_mod_ell'][p] != old_f['ap'][p]:
            return False
    return True

