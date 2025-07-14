/*
  Dual pairs for Galois representations with big image in GL_2(F_l)
  for l = 3 and l = 5
*/

\\ following functions taken from dual-pairs-experimental/etale-algebras.gp

\\ tensor product of two matrices
mattensor(A, B) =
{
   matconcat(apply(a -> a*B, A));
}

multiplication_tensor_split(M) =
{
   my(d = matsize(M)[1]);
   vector(d, i, M[, (i-1)*d+1..i*d]);
}

\\ multiplication tensor of Q[x]/(f) w.r.t. the power basis
multiplication_tensor(f) =
{
   my(d = poldegree(f), m = matcompanion(f));
   concat(powers(m, d - 1, matid(d)));
}

\\ matrix of the algebra endomorphism Q[x]/(f) -> Q[y]/(g)
\\ defined by sending x to a (assuming g | f(a))
algebra_homomorphism_matrix(f, g, a) =
{
   my(d = poldegree(f), e = poldegree(g));
   a = Mod(a, g);
   Mat([Colrev(liftpol(b), e) | b <- powers(a, d - 1)]);
}

\\ end of functions taken from dual-pairs-experimental/etale-algebras.gp

\\ express an element of L = (Q[y]/(f))[x]/(g) on the basis
\\ (1, y, ..., y^(l^2-2), x, y*x, ..., y^(l^2-2)*x^(l^2-l-1))
algtobasis_rel(g, a) =
{
   a = liftpol(Mod(a, g));
   concat([Colrev(polcoef(a, i, 'x), l^2 - 1) | i <- [0..l^2 - l - 1]]);
}

\\ assume h square-free
algisincl(K, h) =
{
   my(F = factor(h)[,1],
      incl = [Mod(nfisincl(K, f), f) | f <- F],
      result = []);
   forvec(j = [[1, #i] | i <- incl],
      result = concat(result, [chinese([incl[k][j[k]] | k <- [1..#F]])]));
   result;
}

\\ construct the algebra B = Q × Q[x]/(f_B) and the pairing Phi
make_dual_algebra_and_pairing(mu) =
{
   T = multiplication_tensor_split(mu~);
   n = #T;
   if(T[1] != matid(n), error("bug: first element is not the identity"));
   found = 0;
   for(i = 2, n,
      gen = T[i];
      chi = charpoly(gen);
      if(issquarefree(chi),
	 found = i;
	 break));
   if(!found,
      error("no algebra generator found"));
   a_0 = nfroots(,chi)[1];
   chi_1 = chi / ('x - a_0);
   [f_B, a_1] = polredabs(chi_1, 1);

   M = matconcat([algebra_homomorphism_matrix(chi, x, a_0);
		  algebra_homomorphism_matrix(chi, f_B, a_1)]);

   pow_gen = vector(n);
   pow_gen[1] = vectorv(n, i, i == 1);
   for(i = 2, n, pow_gen[i] = gen * pow_gen[i - 1]);

   Phi = Mat(pow_gen) * M^-1;
   [f_B, Phi];
}

dual_pair_from_big_image_field(f) =
{
   l = sqrtint(poldegree(f) + 1);
   K = nfinit(subst(f, 'x, 'y));

   /*
     The coordinate ring of the group scheme is A = Q × K.
     The field K = Q[y]/(f) has degree l^2 - 1 and has
     automorphism group of order l - 1.
     The algebra L = (Q[y]/(f))[x]/(g) has degree l^2 - l over K
     and absolute degree (l^2 - 1)(l^2 - l).
   */
   Aut_K = nfgaloisconj(K);
   g = f / vecprod([x - Mod(u, K.pol) | u <- Aut_K]);

   \\ projection maps A → Q, A → K
   proj_0 = matconcat([Mat(1), matrix(1, l^2 - 1)]);
   proj_1 = matconcat([matrix(l^2 - 1, 1), matid(l^2 - 1)]);
   id = matid(l^2);

   \\ inclusion Q → K, multiplication K ⊗ K → K
   unit_K = algebra_homomorphism_matrix(x, f, 0);
   mul_K = multiplication_tensor(f);

   \\ antipode on K
   order_2 = [u | u <- Aut_K, u != y && nfgaloisapply(K, u, u) == y][1];
   iota_K = algebra_homomorphism_matrix(K.pol, K.pol, order_2);
   if(iota_K^2 != matid(l^2 - 1), error("bug: iota_K does not have order 2"));

   \\ "compositum" map K ⊗ K → L
   comp = Mat(concat([[algtobasis_rel(g, y^j * x^i)
		       | j <- [0..l^2 - 2]] | i <- [0..l^2 - 2]]));

   \\ determine the subalgebra Lsym of L fixed under swapping x and y
   \\ TODO: the element x + y might not always generate Lsym
   if(!issquare(norm(subst(charpoly(Mod(x + y, g)), 'x, 'z)), &h),
      error("polynomial not a square"));
   if(!issquarefree(h),
      error("polynomial not square-free: x + y does not generate Lsym"));
   Lsym_to_L = Mat([algtobasis_rel(g, (x + y)^i) |
		    i <- [0..(l^2 - 1)*(l^2 - l)/2 - 1]]);

   \\ inclusions K → L factoring via Lsym
   inclusions = [Lsym_to_L * algebra_homomorphism_matrix(f, h, incl)
		 | incl <- algisincl(K, h)];

   if(l == 3,
      \\ only possibility for the scalar multiplication maps [2]
      scalar_possibilities = [[iota_K]],
      l == 5,
      \\ multiplication maps by 2 and 3 on K, in unknown order
      mult = [algebra_homomorphism_matrix(K.pol, K.pol, u)
	      | u <- Aut_K, u != y && u != order_2];
      \\ possibilities for the scalar multiplication maps [2, 3, 4]
      scalar_possibilities = [[mult[1], mult[2], iota_K],
			      [mult[2], mult[1], iota_K]]);

   foreach(scalar_possibilities, scalars,
      \\ isomorphism A ⊗ A → Q × K^{l+1} × L
      isom = matconcat([mattensor(proj_0, proj_0),
			mattensor(proj_0, proj_1),
			mattensor(proj_1, proj_0),
			mul_K * mattensor(proj_1, proj_1),
			matconcat([mul_K * mattensor(proj_1, s * proj_1)
				   | s <- scalars]~),
			comp * mattensor(proj_1, proj_1)]~);

      \\ Try maps K → Lsym until we find one that gives a valid Hopf
      \\ algebra.  We only need to check coassociativity, all other
      \\ conditions hold by construction.
      foreach(inclusions, incl,
	 M = matconcat([proj_0,
			proj_1,
			proj_1,
			matconcat([s * proj_1 | s <- scalars]~),
			unit_K * proj_0,
			incl * proj_1]~);
	 mu = matsolve(isom, M);
	 if(mattensor(mu, id) * mu == mattensor(id, mu) * mu,
	    [f_B, Phi] = make_dual_algebra_and_pairing(mu);
	    return([[x, f], [x, f_B], Phi]))));
   error("no valid Hopf algebra found");
}

to_lmfdb_format(D) =
{
   [A, B, Phi] = D;
   d = denominator(Phi);
   [apply(Vecrev, A), apply(Vecrev, B), [d, [v~ | v <- d*Phi~]]];
}
