# Mod-l-galois-representations

In this repository we store raw data for the [mod ell Galois representations]  for the [LMFDB](https://github.com/LMFDB/lmfdb).

Each line in a data file will contain a single comma-separated list (enclosed in square brackets) correspond to one irreducible mod-ell Galois representation. The elements of the list are as follows, in this exact order:

* **"numberfield_label"** (string): label of the field whose absolute Galois group is being represented
* **dimension** (integer): the dimension of the representation
* **[finitefield_char, finitefield_degree]** (list of two integers): [ell, degree of rep's ground field over F_ell]
* **conductor** (integer): prime-to-ell Artin conductor of the rep
* **primes_dividing_conductor** (list of ints): [p: p prime dividing conductor]
* **weight** (integer): "weight" of the representation, integer modulo ell - 1 determined by the determinant of the rep
* **absirred_boolean** (int): 0 or 1 depending on whether the rep is absolutely irreducible
* **"type"** (string): "lin" = linear, "orth" = orthogonal, or "sym" = symplectic 
* **"image_type"** (string): "big" = contains SL_n; rest is TBD
* **"image_label"** (string): if known, label of group cf. Sutherland; if not "". 
* **"image_attribute"** (string): "solvable" or TBD
* **image_order** (integer): size of the image if known; if unknown "". 
* **degree_proj_field** (integer): degree of field of definition of the projective representation over F_ell
* **"projective_type"** (string): "A5" or similar
* **"projective_label"** (string): if known
* **bad_prime_list** (list of various attributes): [[int(p), "polynomial", "type", int(order), int(proj_order)] for p dividing conductor]
* **good_prime_list** (list of various attributes): [[int(p), ["tr1", "tr2", ...], "factored_char_poly", int(order), int(proj_order) for p prime < 100 not dividing N*ell]
* **"poly_ker"** (string): polynomial carving out the number field fixed by the kernel of rep
* **"poly_proj_ker"** (string): same, but for the projective rep
* **[related_objects]** (list of pairs of strings): [["object_type", "object_label"] for objects related to rep]

We insist that all strings be enclosed in straight double quotes, "like so"; unknown strings should look like "".

No example currently available.
```

```

This Repository contains the following files:

REST TO BE INSERTED 
