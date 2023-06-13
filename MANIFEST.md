# Mod-l-galois-representations

Where is what?

## code/mfmodell ##

Contains Sage + Magma code to compute mod ell reductions of modular forms

- code/mfmodell/modell.py:  run on legendre.mit.edu in lmfdb directory
			  see comments in code for detailed instructions

- code/mfmodell/bigspaces.m: requires mf.m from CMFs/magma.  To run this, add a symlink to it in CMFs/magma/ and run magma there.

## data/mfmodell ##

Contains data output from modell.py: data/mfmodell/<ell> for ell=0,2,3,5,7 contains data in characteristic ell

ell=2

- data/mfmodell/2/1-1000.txt contains 12739 mod 2 forms from levels 1-1000 (1960 excluding repeats)
- data/mfmodell/2/x1-1000.txt contains 1368 undecomposed spaces of dim>50, not processed for mod 2 forms
- data/mfmodell/2/mod2bb1000.txt contains BlackBox output for 9609 mod 2 reductions (includes repeats)
    one line each saying "reducible" or "S3" or "C3" with cubics defining splitting fields, e.g.

    5.4.a.a:reducible
    11.2.a.a:S3:x^3-x^2-x-1
    49.4.a.c:C3:x^3-7*x-7

ell=3

- data/mfmodell/3/1-1000.txt contains 14293 mod 3 forms from levels 1-1000 (8008 excluding repeats)
- data/mfmodell/3/x1-1000.txt contains  2254 undecomposed spaces of dim>50, not processed for mod 3 forms
- data/mfmodell/3/mod3bb1000nr.txt contains BlackBox output for 1550 mod 3 reductions (excludes repeats)
   one line each saying "reducible" or linear and projective image if
   irreducible with polynomials defining splitting fields, e.g.

    5.4.a.a:GL(2,3):x^8-4*x^7+7*x^6-7*x^5+4*x^4-x^3-4*x^2+4*x-1:S4:x^4-x^3+4*x-1
    7.3.b.a:Ns:(x^4-x^3-3*x^2-x+1)*(x^4-x^3+2*x+1):V4:(x^2-21)*(x^2+3)
    7.4.a.a:GL(2,3):x^8-6*x^4-x^2-3:S4:x^4-x^3-3*x^2-7*x+1
    8.3.d.a:reducible

ell=5

- data/mfmodell/5/1-1000.txt contains 24851 mod 5 forms from levels 1-1000 (21177 excluding repeats)
- data/mfmodell/5/x1-1000.txt contains  3145 undecomposed spaces of dim>50, not processed for mod 5 forms

ell=7

- data/mfmodell/7/1-1000.txt  contains 32761 mod 7 forms from levels 1-1000 (30425 excluding repeats)
- data/mfmodell/7/x1-1000.txt  contains 4213 undecomposed spaces of dim>50, not processed for mod 7 forms

ell=0

- Contains magma output for many spaces of dimention 21-50 used in
  computing the above, but I stopped adding these to the repo as they
  were getting too big (>5G in all).
