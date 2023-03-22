# Mod-l-galois-representations

Description of columns for the LMFDB table [modlgal_reps](https://beta.lmfdb.org/api/modlgal_reps).

Currently we are only considering representations of Gal_Q whose base ring is the field with l elements, but the schema is more general.

* **numberfield_label** (text): label of the field whose absolute Galois group is represented (Q=1.1.1.1);
* **dimension** (integer): the dimension of the representation;
* **algebraic_group** (text): the algebraic group for the codomain (e.g. "GL", "GSp");
* **base_ring_characteristic** (integer): characteristic of the base ring (the value of l);
* **base_ring_order** (integer): order of the base ring;
* **base_ring_is_field** (boolean): true iff base ring is a field;
* **conductor** (integer): prime-to-ell Artin conductor of the representation;
* **conductor_primes** (integer[]): list of primes dividing the conductor;
* **conductor_squarefree** (boolean): true iff the conductor is squarefree;
* **determinant_label** (text): label of the mod ell Galois character which is the determinant of the representation;
* **cyclotomic_exponent** (integer): exponent of the cyclotomic character of the representation in the determinant of the representation, an integer modulo (l - 1);
* **image_label** (text): label of image (identitifies entry in gl1zhat, gl2zhat, gsp4zhat, tables);
* **image_type** (text): human sensible desscription of the image, e.g "big", "cyclic" etc;
* **image_index** (integer): index of the image;
* **image_order** (integer): order of the image;
* **is_solvable** (boolean): true iff image is a solvable group
* **is_surjective** (boolean): true iff image equals codomain
* **is_irreducible** (boolean): true iff image is irreducible (currently always true, only set when base ring is a field)
* **is_absolutely_irreducible** (boolean): true if image is absolutely irreducible (only set when base ring is a field)
* **projective_label** (string): TBD;
* **projective_type** (string): human sensible descriptor of the projective image, e.g. "A5" or "trivial" similar;
* **projective_is_surjective** (boolean): true iff projective rep is surjective;
* **good_prime_list** (jsonb): list of lists of the form 

      [p, trace, determinant, charpoly, factorization_charpoly, order, projective_order, representative_matrix]
   
   where 
   - p is int(p) or string if p is a prime in a number field different from Q;
   - trace, determinant are lists of integers which represent the corresponding elements in the coefficient field;
   - charpoly (list of lists): coefficients of the characteristic polynomial of the Frobenius matrix at p represented as above;
   - factorization_charpoly (list of pairs): every pair is a lists consisting of a list of coefficients of an irreducible factor of charpoly and its exponent in the factorization;
   - order, projective order (integers): order, projective order of the Frobenius matrix at p;
   - representative_matrix (list of lists): list of length dimension^2 consisting of coefficients of a representative for the conjugacy class (inside the image) of a Frobenius matrix at p;

for primes p< 100 (we will decide later about number fields) **unramified** (including finitefield_char is the representation is unramified at finitefield_char);
* **real_place_list** (jsonb): same as for good_prime_list where the first entry is a 0-based index into the ordered list of real roots of the defining polynomial of numberfield_label;         
* **bad_prime_list** (jsonb): list of lists. For each ramified prime give a list of 4 elements:

      [p, inertia_invariants, inertia_coinvariants, type] 
 
 where p is a prime, inertia_invariants and inertia_coinvariants are lists of the form
 
          [trace, determinant, charpoly, factorization_charpoly, order, projective_order, representative_matrix]

 everything represented with the same format as good_prime_list, and type is a list of descriptors. For 2-dimensional representation type = [is_reducible, is_decomposable] (list of booleans).
* **kernel_polynomial** (integer[]): list of integer coefficients of polynomial defining canonical sibling (chosen as with Artin stem fields) whose splitting field is the fixed field of the kernel of the representation;
* **projective_kernel_polynomial** (integer[]): kernel polynomial for the projective representation
* **related_objects** (text[]): each pair is of the form ["object_type", "object_label"] for objects related to representation.
