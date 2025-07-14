/*
  Abelian (one-dimensional) representations of the absolute
  Galois group of the rational numbers over finite fields

  Simplified version (only for prime fields) of the code at
  https://gitlab.com/pbruin/abgalrep
  Returns the dual pair in LMFDB format.
*/

\\ Return a pair [d, e], where d is the gcd of the elements of v
\\ and e is a column vector such that v*e = d.
vecgcdext(v) =
{
   my(H, U);
   if(v == [],
      return([0, []]));
   [H, U] = mathnf(v, 1);
   [H[1, 1], U[,length(v)]];
}

\\ Return g in (Z/mZ)^* mapping to 1/o (= "zeta") under chi.
character_split(D, chi) =
{
   my(o = charorder(D, chi),
      ord = D[4][5],  \\ orders of Conrey generators
      r = length(ord),
      c = if(type(chi) == "t_COL", chi, znconreylog(D, chi)),
      q = vector(r, i, c[i]/ord[i]),
      d, e, g);
   if(denominator(q) != o,
      error("denominator != o"));
   [d, e] = vecgcdext(concat(q * o, [o]));
   g = znconreyexp(D, e[1..r]);
   if(chareval(D, chi, g) != (1 % o)/o,
      error("wrong generator found"));
   g;
}

good_prime(m, g) =
{
   forprime(p = 2,,
      if(p % m == g, return(p)));
}

\\ Given a subgroup H of (Z/mZ)^* (as a list of generators),
\\ return [K, alpha_g(x)], where K is the fixed field of H
\\ and alpha_g(x) is the image of x under the automorphism
\\ of Q(zeta_m) sending zeta_m to zeta_m^g.
galoissubcyclo_with_aut(D, H, g, v='x) =
{
   my(m = D.mod,
      f = polredabs(galoissubcyclo(m, H, 0, v)),
      K = nfinit(f),
      G = galoisinit(K),
      p = good_prime(m, g),
      s = idealfrobenius(K, G, idealprimedec(K, p)[1]),
      h = galoispermtopol(G, s));
   [K, Mod(h, f)];
}

\\ Return the orbit of x under the automorphism
\\ aut of L (ordered cyclically).
nfgaloisorbit(L, aut, x) =
{
   my(orbit = [x]);
   while((x = nfgaloisapply(L, aut, x)) != orbit[1],
         orbit = concat(orbit, [x]));
   orbit;
}

\\ convert H from HNF format to a list of generators
subgroup_generators(D, H) =
{
   my(g = Mod(D.gen, D.mod));
   [prod(i = 1, #g, g[i]^v[i]) | v <- H];
}

\\ Return the character of Gal(L_chi/Q) corresponding to
\\ the Dirichlet character chi with values in <z>.
galois_character_as_table(D, chi, z, z_order, var='x) =
{
   my(H = subgroup_generators(D, charker(D, chi)),
      \\ directly doing H = charker(D, chi) may give wrong results
      \\ when fed to galoissubcyclo, since e.g. znstar(55) and
      \\ znstar(55, 1) have different SNF generators!
      o = charorder(D, chi),
      g = character_split(D, chi),  \\ element of (Z/mZ)^* mapping to z
      L, aut, orb);
   [L, aut] = galoissubcyclo_with_aut(D, H, g, var);
   \\ aut is the automorphism corresponding to z (and g) via chi
   orb = nfgaloisorbit(L, aut, Mod(variable(L), L.pol));
   [L, Mat([orb~, powers(z^(z_order/o), o - 1)~])];
}

\\ l-cyclotomic character mod l
cyclotomic_character(Dl) =
{
   if(Dl.no == 1, 1, znconreyexp(Dl, [1]));
}

\\ Return a list of all orbits of z in F_l^* acting on F_l,
\\ ordered cyclically: [x, z*x, z^2*x, ...]
zeta_orbits(l, z) =
{
   my(orbits = [], X, x, S);
   X = Set(Mod([0..l-1], l));
   while(X != [],
         x = X[1];
         S = [x];
         while((x *= z) != S[1],
               S = concat(S, [x]));
         orbits = concat(orbits, [S]);
         X = setminus(X, Set(S)));
   orbits;
}

\\ incl = root of K.pol in L
\\ return the matrix of the inclusion map
incl_matrix(K, L, incl) =
{
   my(d = poldegree(K.pol),
      e = poldegree(L.pol));
   incl = Mod(incl, L.pol);
   Mat([Colrev(x, e) | x <- liftpol(powers(incl, d - 1))]);
}

\\ incl = root of K.pol in L, or matrix of inclusion map
\\ apply the corresponding inclusion map to x
apply_incl(K, L, incl, x) =
{
   x = liftpol(x);
   if(type(incl) == "t_MAT",
      Mod(Polrev(incl * Colrev(x, poldegree(K.pol)), variable(L.pol)), L.pol),
      subst(x, variable(K), Mod(incl, L.pol)));
}

char_dual_pair_2(D1, chi1, D2, chi2, l, z, z_order) =
{
   my(order1 = charorder(D1, chi1),
      order2 = charorder(D2, chi2),
      orbits1 = zeta_orbits(l, z^(z_order/order1)),
      orbits2 = zeta_orbits(l, z^(z_order/order2)),
      no1 = #orbits1, no2 = #orbits2,
      L1, table1, L2, table2, Mdata, f, incl1, incl2,
      zeta_l, zeta_powers, P0, P, Q0, Q, W, e,
      polys1, polys2, Phi, d);
   [L1, table1] = galois_character_as_table(D1, chi1, z, z_order, 'x);
   [L2, table2] = galois_character_as_table(D2, chi2, z, z_order, 'y);
   Mdata = polcompositum(subst(L1.pol, 'x, 'z),
                         subst(L2.pol, 'y, 'z), 1)[1];
   f = Mdata[1];
   incl1 = incl_matrix(L1, f, Mdata[2]);
   incl2 = incl_matrix(L2, f, Mdata[3]);
   zeta_l = Mod(nfisincl(polcyclo(l), f)[1], f);
   zeta_powers = powers(zeta_l, l - 1);
   \\ matrix of the isomorphism A_{chi1,S,M} -> M^S for S non-trivial
   P0 = Mat([[apply_incl(L1, f, incl1, nfgaloisapply(L1, a, b))
              | b <- powers('x, poldegree(L1.pol) - 1)]
             | a <- table1[,1]]~);
   \\ total matrix
   P = matconcat(matdiagonal(vector(no1, k,
                                    if(#orbits1[k] == 1, Mat(1), P0))));
   \\ matrix of the isomorphism A_{chi2,S,M} -> M^S for S non-trivial
   Q0 = Mat([[apply_incl(L2, f, incl2, nfgaloisapply(L2, a, b))
              | b <- powers('y, poldegree(L2.pol) - 1)]
             | a <- table2[,1]]~);
   \\ total matrix
   Q = matconcat(matdiagonal(vector(no2, k,
                                    if(#orbits2[k] == 1, Mat(1), Q0))));
   orbits1 = concat(orbits1);
   orbits2 = concat(orbits2);
   \\ matrix of the pairing A \times B -> Q
   W = matrix(#orbits1, #orbits2, i, j,
              e = orbits1[i] * orbits2[j];
              zeta_powers[liftint(e) + 1]);
   polys1 = concat(['x], vector(no1 - 1, i, L1.pol));
   polys2 = concat(['y], vector(no2 - 1, j, L2.pol));
   Phi = liftpol(P~ * W * Q) / #orbits1;
   d = denominator(Phi);
   [apply(Vecrev, polys1), apply(Vecrev, polys2),
    [d, [v~ | v <- d * Phi~]]];
}

\\ Dual pair of algebras attached to a Dirichlet character over F_l
\\ (where l does not divide the order of chi)
\\ z must be a root of unity such that its order z_order is
\\ divisible by lcm(order(chi), l - 1), and such that
\\ z^(z_order/(l - 1)) equals the standard primitive root in Z/lZ.
char_dual_pair_1(D, chi, l, z, z_order) =
{
   my(Dl = znstar(l, 1),
      D2 = znstar(lcm(D.mod, l), 1),
      chi2 = chardiv(D2, zncharinduce(Dl, cyclotomic_character(Dl), D2),
                     zncharinduce(D, chi, D2)));
   \\ [D, chi] = znchartoprimitive(D, chi);
   \\ [D2, chi2] = znchartoprimitive(D2, chi2);
   char_dual_pair_2(D, chi, D2, chi2, l, z, z_order);
}

abgalrep(n, chi, l) =
{
   my(D = znstar(n, 1),
      o = charorder(D, chi),
      d = znorder(Mod(l, o)),
      z = znprimroot(l));
   char_dual_pair_1(D, chi, l, z, l^d - 1);
}
