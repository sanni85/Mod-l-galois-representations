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
  format: one per line, fields separated by :, listed in code/mfmodell/modell.py's nf_to_string().
  NB if a newform has several different mod-ell reductions this is not recorded in the labelling in this file

- data/mfmodell/2/x1-1000.txt contains 1368 undecomposed spaces of dim>50, not processed for mod 2 forms
- data/mfmodell/2/mod2bb1000.txt contains BlackBox output for 12739 mod 2 reductions:

    - one line each with fields separated by ":".
    - First field is the modular form label with additional "-1", "-2",
      etc to distinguish multiple reductions of the same modular form.
    - Second field is "reducible" or "S3" or "C3".
    - Third field (for irreducibles) is a cubic defining the splitting field.

    e.g.

    5.4.a.a-1:reducible
    11.2.a.a-1:S3:x^3-x^2-x-1
    49.4.a.c-1:C3:x^3-7*x-7

- data/mfmodell/2/mod2dim2reps.yxy is a dump from the database modlgal table of rows with ell=2 and dim=2.
- data/mfmodell/2/label_map_2.txt has one line for each database entry
  for ell=2 showing which modular forms match it (a_p for good p<100
  and splitting fields match).
- data/mfmodell/2/reverse_label_map_2.txt: inverse of previous.


ell=3

- data/mfmodell/3/1-1000.txt contains 14293 mod 3 forms from levels 1-1000 (8008 excluding repeats)
- data/mfmodell/3/x1-1000.txt contains  2254 undecomposed spaces of dim>50, not processed for mod 3 forms
- data/mfmodell/3/mod3bb1000.txt contains BlackBox output for 24569 mod 3 reductions:

    - one line each with fields separated by ":".
    - First field is the modular form label with additional "-1", "-2",
      etc to distinguish multiple reductions of the same modular form.
    - 2nd field is "reducible" or "GL(2,3)" or "Nn" or "Ns" (image).
    - 3rd field (for irreducibles) is an octic defining the splitting field.
    - 4th field (for irreducibles) is "S4" or "V4" or "D4" (projective image)
    - 5th field (for irreducibles) is a quartic defining the projective splitting field.

    5.4.a.a-1:GL(2,3):x^8-4*x^7+7*x^6-7*x^5+4*x^4-x^3-4*x^2+4*x-1:S4:x^4-x^3+4*x-1
    7.3.b.a-1:Ns:(x^4-x^3-3*x^2-x+1)*(x^4-x^3+2*x+1):V4:(x^2-21)*(x^2+3)
    7.4.a.a-1:GL(2,3):x^8-x^7-3*x^6+2*x^5+4*x^4+3*x^3-5*x^2-7*x-3:S4:x^4-x^3-3*x^2-7*x+1
    8.3.d.a-1:reducible

- data/mfmodell/3/mod3dim2reps.yxy is a dump from the database modlgal table of rows with ell=3 and dim=2.
- data/mfmodell/3/label_map_3.txt has one line for each database entry
  for ell=3 showing which modular forms match it (a_p for good p<100
  and splitting fields and projective splitting fields match).
- data/mfmodell/3/reverse_label_map_3.txt: inverse of previous.

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
