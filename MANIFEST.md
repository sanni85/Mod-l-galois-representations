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

- data/mfmodell/2/1-1000.txt contains  9609 mod 2 forms from levels 1-1000 (1550 excluding repeats)
- data/mfmodell/2/x1-1000.txt contains  2859 undecomposed spaces, not processed for mod 2 forms
- data/mfmodell/2/mod2bb1000.txt contains BlackBox output for 9609 mod 2 reductions (includes repeats)
    one line each saying "reducible" or "S3" or "C3" with cubics defining splitting fields
- mod2bb1000_concise.txt concise version of previous

ell=3

- data/mfmodell/3/1-1000.txt contains 11120 mod 3 forms from levels 1-1000 (6687 excluding repeats)
- data/mfmodell/3/x1-1000.txt contains  4512 undecomposed spaces, not processed for mod 3 forms
- data/mfmodell/3/mod3bb1000nr.txt contains BlackBox output for 1550 mod 3 reductions (excludes repeats)
    one line for "reducible" or two lines for linea and projective image if irreducible with polynomials defining splitting fields
- mod3bb1000nr_concise.txt concise version of previous

ell=5


- data/mfmodell/5/1-1000.txt contains 19528 mod 5 forms from levels 1-1000 (6758 excluding repeats)
- data/mfmodell/5/x1-1000.txt contains  6758 undecomposed spaces, not processed for mod 5 forms

ell=7

- data/mfmodell/7/1-1000.txt  contains 25890 mod 7 forms from levels 1-1000 (24419 excluding repeats)
- data/mfmodell/7/x1-1000.txt  contains 9811 undecomposed spaces, not processed for mod 7 forms


