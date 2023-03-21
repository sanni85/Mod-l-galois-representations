# Mod-l-galois-representations

Where is what?

code/mfmodell contains Sage + Magma code to compute mod ell reductions of modular forms

code/mfmodell/modell.py:  run on legendre.mit.edu in lmfdb directory
			  see comments in code for detailed instructions

code/mfmodell/bigspaces.m: requires mf.m -- where is that?

data/mfmodell contains data output from modell.py
data/mfmodell/<ell> for ell=0,2,3,5,7 contains data in characteristic ell

data/mfmodell/2/1-1000.txt contains  9609 mod 2 forms from levels 1-1000 (1550 excluding repeats)
data/mfmodell/2/x1-1000.txt contains  2859 undecomposed spaces, not processed for mod 2 forms
data/mfmodell/2/mod2bb1000.txt contains BlackBox output for 9609 mod 2 reductions (includes repeats)
    one line each saying "reducible" or "S3" or "C3" with cubics defining splitting fields

data/mfmodell/3/1-1000.txt contains 11120 mod 3 forms from levels 1-1000 (6687 excluding repeats)
data/mfmodell/3/x1-1000.txt contains  4512 undecomposed spaces, not processed for mod 3 forms
data/mfmodell/3/mod3bb1000nr.txt contains BlackBox output for 1550 mod 3 reductions (excludes repeats)
    one line for "reducible" or two lines for linea and projective image if irreducible with polynomials defining splitting fields

data/mfmodell/5/1-1000.txt contains 19528 mod 5 forms from levels 1-1000 (6758 excluding repeats)
data/mfmodell/5/x1-1000.txt contains  6758 undecomposed spaces, not processed for mod 5 forms

data/mfmodell/7/1-200.txt  contains  5573 mod 7 forms from levels 1-1000 (5251 excluding repeats)
data/mfmodell/7/x1-200.txt  contains  1051 undecomposed spaces, not processed for mod 7 forms


