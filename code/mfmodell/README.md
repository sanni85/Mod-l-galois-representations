# Code for computing mod-ell reductions of modular forms #

## modell.py ##

For details see comments in modell.py.  It is a Sage script which contains functions

- run(levels, ells, max_dim, verbose)
- get_forms(N,k,ell)
- get_form_data_from_file()
- data_output()
- extra_output()

as well as various utilities / subroutines.

run() computes all mod ell reductions (for each ell in ells) of all
newforms of level N (for each N in levels) and weight k (for 2 <= k <=
max(ell+1,4)) by calling get_forms(N,k,ell) for each triple.

get_forms() accesses the LMFDB tables mf_newforms and mf_hecke_nf, and
returns two lists, one of newforms with Hecke field data and at least
one mod ell reduction, the other a list of labels of newforms with no
Hecke field information in the database. If the second list is
nonempty after reading from the database, run() attempts to read the
missing newform data from a data file (with filename = newform label),
using get_form_data_from_file(), but only if the newform dimension is
<= max_dim.  This is the only place where the parameter max_dim is
used).

To create a suitable data file, use the Magma script bigspaces.m.

Rationale: the LMFDB only contains Hecke field data for newforms of
dimension at most 20.  In early runs (with a small upper bound on the
level) this left only a few gaps, which we wanted to fill by running
those spaces separately.  Now that we have covered all levels up to
1000 this is not practical.  For example, running the 5 spaces
95.2.p.a, 97.4.c.a, 97.4.d.a, 99.4.j.a, 88.6.c.a of dimension 48 or 50
takes a few hours.

The rest of the top-level function run() just accumulates all the mod
ell reductions into a dict which is returned.  This includes, for each
mod ell reduction, an index i, meaning that the same list of (ap mod
ell) has been seen i-1 times before.  Thus to obtain a list with no
repeats we can select those with i=1.

To produce the output files use the functions data_output() and
extra_output() with the dict returned by run() as first parameter and
a filename as second parameter.

## bigspaces.m ##

This is just a front end to magma code in the CMFs repository, from
which it loads mf.m.  I run this by putting a copy (or symlink) into
CMFs/magma.   It contains two procedures:

- OneSpace(label: Coeffs:=100, DegreeBound:=0, Detail:=0)
- DoAll(fname :  spaces_done := [], DegreeBound:=0, Detail:=0)

OneSpace() takes a newform label with 3 or 4 components and calls
NewSpaceData() to compute all the required data for the newform,
outputting to a file, suitable for reading by the Sage function
get_form_data_from_file().

DoAll() does the same for all newforms with labels in the file given
as first parameter.  If more than one newform in the file are in the
same space (i.e. same N,k,c components of the label N.k.c.i) then
OneSpace() is only called the first time.
