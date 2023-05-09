#!/usr/local/bin/sage -python
# -*- coding: utf-8 -*-
r""" Make 1-dim mod l representations

"""

import re
import os
import sys

HOME=os.path.expanduser("~")
sys.path.append(os.path.join(HOME, 'lmfdb'))

from lmfdb import db
from sage.all import ZZ, primitive_root, prime_range, GF, DirichletGroup, QQ, PolynomialRing

reps = db.modlgal_reps
chars = db.char_dir_values

# precompute cyclopowers
cyclopowers={'2':[0,0],'3':[0,0,1], '5':[0,0,1,3,2], '7':[0,0,2,1,4,5,3]}

Qx=PolynomialRing(QQ, 'x')
# Label will be ell.1.N.c where 1=dim, N=conductor, c=Conrey #

# Input is character chi, and the ell we are reducing
def buildrep(chi, ell):
  conrey = max(chi.conrey_number(), 1)
  N = chi.conductor()
  charlabel = rf'{N}.{conrey}'
  charlabel1 = rf'{N}-{conrey}'
  N1 = ZZ(N/ell**N.valuation(ell))
  F = GF(ell)
  primroot = primitive_root(ell)

  #charlabel = ent['label']
  cycloexp = discrete_log(F(conrey), F(primroot), ell-1) if ZZ(ell).divides(N) else 0
  badprimes = N1.prime_divisors()
  newlabel = rf'{ell}.1.{N1}.{charlabel1}'
  order = chi.order()
  # Need kernel field
  kerpol = chi.fixed_field_polynomial()
  kerpol = Qx(pari(kerpol).polredabs()).coefficients(sparse=False)
  kerpol = [ZZ(z) for z in kerpol]

  plist = [p for p in prime_range(100) if not p.divides(ell*N)]
  goodplist=[[p, [ZZ(chi(p))]] for p in plist]
  mygen = 0
  for j in range(len(plist)):
    if F(goodplist[j][1][0]).multiplicative_order()==order:
        mygen = plist[j]
        break
  assert mygen > 0
  rec = {'algebraic_group': 'GL', 
         'bad_prime_list': [[z,''] for z in badprimes],
         'base_ring_characteristic': ell,
         'base_ring_is_field': True, 
         'base_ring_order': ell, 
         'conductor': N1,
         'conductor_primes': badprimes, 
         'conductor_is_squarefree': N1.is_squarefree(),
         'conductor_num_primes': len(badprimes),
         'cyclotomic_exponent': cycloexp, 
         'determinant_label': newlabel, 
         'dimension': 1,
         'good_prime_list': goodplist,
         'image_index':ZZ((ell-1)/order), 
         'image_label': '$C_{%d}$'%order, 
         'image_order': order, 
         'image_type': 'big' if order> 1 else 'trivial', 
         'is_absolutely_irreducible': True,
         'is_irreducible': True, 
         'is_solvable': True, 
         'is_surjective': order == ell-1, 
         'kernel_polynomial': kerpol,
         'label': newlabel, 
         'projective_is_surjective': ell==2, 
         'projective_kernel_polynomial':[0,1], 
         'projective_type': 'big' if order> 1 else 'trivial',
         'top_slope_rational': '0' if N.valuation(ell)==0 else '1', 
         'top_slope_real':0 if N.valuation(ell)==0 else 1, 
         'generating_primes': [mygen], 
         'related_objects': [['Dirichlet', charlabel]],
         'image_abstract_group': rf'{order}.1',
         'projective_image_abstract_group': '1.1',
         'frobenius_matrices': [z[1] for z in goodplist]}
  return rec

def allc(c, ell):
    result = []
    if not ZZ(ell).divides(c):
      for N in [c, c*ell]:
        F = GF(ell)
        primroot = primitive_root(ell)
        #dets = [Matrix(F,n,z).det() for z in frobs]
        DG=DirichletGroup(N, F, zeta=primroot)
        call = [z for z in DG if z.order().divides(ell-1)]
        call = [z for z in call if z.is_primitive()]
        result += [buildrep(z,ell) for z in call]
    return result

def flatten(l):
    return [item for sublist in l for item in sublist]

#outrecs.append(data)
#reps.insert_many(outrecs)

