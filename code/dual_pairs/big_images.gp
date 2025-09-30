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
algtobasis_rel(g, a, l) =
{
   a = liftpol(Mod(a, g));
   concat([Colrev(polcoef(a, i, 'x), l^2 - 1) | i <- [0..l^2 - l - 1]]);
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

   M = matconcat([algebra_homomorphism_matrix(chi, 'x, a_0);
		  algebra_homomorphism_matrix(chi, f_B, a_1)]);

   pow_gen = vector(n);
   pow_gen[1] = vectorv(n, i, i == 1);
   for(i = 2, n, pow_gen[i] = gen * pow_gen[i - 1]);

   Phi = Mat(pow_gen) * M^-1;
   [['x, f_B], Phi];
}

find(v, P) = for(j = 1, #v, if(P(v[j]), return(j)));

candidate_incl(K, Aut_K, p) =
{
   my(L = rnfequation(K, p, 1),
      \\ compute all inclusions K -> L = K[x]/(p)
      \\ (via the corresponding absolute number field)
      incl_abs = nfisincl(K, L[1]),
      incl_rel = subst(incl_abs, 'x, Mod('x, p) + L[3] * Mod('y, K.pol)),
      \\ find the ones that are do not correspond to just taking of the
      \\ coordinates multiplied by a scalar
      triv = concat(Aut_K, subst(Aut_K, 'y, 'x)),
      candidates = [i | i <- incl_rel, !find(triv, t -> i == t)]);

   \\ find candidate addition laws such that the maps (P, Q) -> (P + Q, Q)
   \\ and (P, Q) -> (P, P + Q) respect L_1
   [i | i <- candidates,
    Mod(substvec(liftpol(p), ['x, 'y], [i, 'y]), p) == 0 &&
    Mod(substvec(liftpol(p), ['x, 'y], ['x, i]), p) == 0];
}

complete_incl(K, i_1, G, scalars) =
{
   my(l = #scalars + 1,
      scalars2 = [subst(u, 'y, 'x) | u <- scalars],
      i = vector(l - 1),
      j, s, t);
   i[1] = i_1;
   for(d = 2, l - 1,
      s = Mod(scalars2[d], G[1]);
      j = find(G, p -> subst(p, 'x, s) == 0);
      t = substvec(liftpol(i[d - 1]), ['x, 'y],
		   Mod(Mod([scalars2[d - 1], i_1], K.pol), G[1]));
      i[d] = substvec(liftpol(t), ['x, 'y],
		      Mod(Mod([scalars2[1/d % l], 'y], K.pol), G[j])));
   iferr(chinese(i), e, 0);
}

good_inclusions(K, g, G, incl_1, scalars) =
{
   my(incl = [complete_incl(K, i_1, G, scalars) | i_1 <- incl_1]);
   \\ find the symmetric ones
   [i | i <- incl, i != 0 &&
    i == Mod(Mod(substvec(liftpol(i), ['x, 'y], ['y, 'x]), g), K.pol)];
}

hopf_algebra_from_big_image_field(f) =
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
   g = f / vecprod(['x - Mod(u, K.pol) | u <- Aut_K]);

   order_2 = [u | u <- Aut_K, u != 'y && nfgaloisapply(K, u, u) == 'y][1];
   if(l == 3,
      scalar_pol = [['y, order_2]],
      l == 5,
      order_4 = [u | u <- Aut_K, u != 'y && u != order_2];
      scalar_pol = [['y, order_4[1], order_4[2], order_2],
		    ['y, order_4[2], order_4[1], order_2]]);

   \\ projection maps A → Q, A → K
   proj_0 = matconcat([Mat(1), matrix(1, l^2 - 1)]);
   proj_1 = matconcat([matrix(l^2 - 1, 1), matid(l^2 - 1)]);
   id = matid(l^2);

   \\ inclusion Q → K, multiplication K ⊗ K → K
   unit_K = algebra_homomorphism_matrix('x, f, 0);
   mul_K = multiplication_tensor(f);

   \\ "compositum" map K ⊗ K → L
   comp = Mat(concat([[algtobasis_rel(g, 'y^j * 'x^i, l)
		       | j <- [0..l^2 - 2]] | i <- [0..l^2 - 2]]));

   G = nffactor(K, g)[,1];

   if(#G == 1,
      \\ determine the subalgebra Lsym of L fixed under swapping x and y
      \\ TODO: the element x + y might not always generate Lsym
      if(!issquare(norm(subst(charpoly(Mod('x + 'y, g)), 'x, 'z)), &h),
	 error("polynomial not a square"));
      if(!issquarefree(h),
	 error("polynomial not square-free: x + y does not generate Lsym"));
      Lsym_to_L = Mat([algtobasis_rel(g, ('x + 'y)^i, l) |
		       i <- [0..(l^2 - 1)*(l^2 - l)/2 - 1]]);
      \\ all inclusions K → L factoring via Lsym
      inclusions = [Lsym_to_L * algebra_homomorphism_matrix(f, h, incl)
		    | incl <- nfisincl(K, h)],
      \\ #G > 1
      incl_1 = candidate_incl(K, Aut_K, G[1]));

   foreach(scalar_pol, scalars,
      scalar_mat = [algebra_homomorphism_matrix(K.pol, K.pol, u)
		    | u <- scalars[2..l-1]];
      \\ isomorphism A ⊗ A → Q × K^{l+1} × L
      isom = matconcat([mattensor(proj_0, proj_0),
			mattensor(proj_0, proj_1),
			mattensor(proj_1, proj_0),
			mul_K * mattensor(proj_1, proj_1),
			matconcat([mul_K * mattensor(proj_1, s * proj_1)
				   | s <- scalar_mat]~),
			comp * mattensor(proj_1, proj_1)]~);

      if(#G > 1,
	 candidates = good_inclusions(K, g, G, incl_1, scalars);
	 inclusions = [Mat([algtobasis_rel(g, incl^i, l) |
			    i <- [0..l^2 - 2]])
		       | incl <- candidates]);

      \\ Try maps K → L until we find one giving a valid Hopf algebra.
      \\ We only need to check coassociativity, all other conditions
      \\ hold by construction.
      foreach(inclusions, incl,
	 M = matconcat([proj_0,
			proj_1,
			proj_1,
			matconcat([s * proj_1 | s <- scalar_mat]~),
			unit_K * proj_0,
			incl * proj_1]~);
	 mu = matsolve(isom, M);
	 if(mattensor(mu, id) * mu == mattensor(id, mu) * mu,
	    return([['x, f], mu]))));
}

dual_pair_from_big_image_field(f) =
{
   my(F, mu, G, Phi);
   [F, mu] = hopf_algebra_from_big_image_field(f);
   [G, Phi] = make_dual_algebra_and_pairing(mu);
   [F, G, Phi];
}

to_lmfdb_format(D) =
{
   [A, B, Phi] = D;
   d = denominator(Phi);
   [apply(Vecrev, A), apply(Vecrev, B), [d, [v~ | v <- d*Phi~]]];
}
