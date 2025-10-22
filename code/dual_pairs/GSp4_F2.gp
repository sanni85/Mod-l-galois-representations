/*
  Dual pairs of algebras for symplectic 4-dimensional irreducible
  representations of Gal(Qbar/Q) over F_2
*/

\r util.gp

default("realprecision", 100);

\\ Reference: D.P. Roberts, Twin sextic algebras, Section 3

duads = concat([[[i, j] | j <- [i+1..6]] | i <- [1..6]]);

synthemes =
{
   [vecextract(duads, i) |
    i <- [[1, 10, 15], [1, 11, 14], [1, 12, 13], [2, 7, 15], [2, 8, 14],
	  [2, 9, 13], [3, 6, 15], [3, 8, 12], [3, 9, 11], [4, 6, 14],
	  [4, 7, 12], [4, 9, 10], [5, 6, 13], [5, 7, 11], [5, 8, 10]]];
}

totals =
{
   [vecextract(synthemes, i) |
    i <- [[1, 5, 9, 11, 13], [1, 6, 8, 10, 14], [2, 4, 8, 12, 13],
	  [2, 6, 7, 11, 15], [3, 4, 9, 10, 15], [3, 5, 7, 12, 14]]];
}

p(r, s) =
{
   my(a = r[s[1][1]], b = r[s[1][2]],
      c = r[s[2][1]], d = r[s[2][2]],
      e = r[s[3][1]], f = r[s[3][2]]);
   a^2 * b^2 * (c * d + e * f) +
   c^2 * d^2 * (e * f + a * b) +
   e^2 * f^2 * (a * b + c * d);
}

y(r, t) = vecsum([p(r, s) | s <- t]);

\\ polynomial defining the twin sextic algebra
\\ TODO: adjust precision if needed
twin(f) =
{
   my(r = polroots(f));
   if(poldegree(f) == 5,
      r = concat([0], r));
   polredabs(bestappr(vecprod(['x - y(r, t) | t <- totals])));
}

\\ subsets of size 0 or 2 of [1, 2, 3, 4, 5, 6]
S = concat([[]], duads);
W = Mat([[(-1)^#setintersect(P, Q) | P <- S]~ | Q <- S]);
Winv = 1/16 * W;  \\ = W^-1

\\ polynomial defining the algebra of unordered pairs of roots of f
pairs_algebra(f) =
{
   my(g = f / ('x - Mod('y, subst(f, 'x, 'y))), k, q, H, h);
   for(k = 0, 2,
      q = 'x + 'y + k * 'x * 'y;
      H = norm(subst(charpoly(Mod(q, g)), 'x, 'z));
      if(!issquare(H, &h),
	 error("not a square"));
      if(issquarefree(h),
	 return([h, q])));
   error("no generator found");
}

\\ TODO: adjust precision if needed
dual_pair_GSp4_F2(f) =
{
   my(R = polroots(f), h1, h, q, T, M, red, F, C, D, Phi);
   [h1, q] = pairs_algebra(f);
   if(poldegree(f) == 5,
      R = concat([0], R);
      h1 = subst(f, 'x, 'z) * h1);
   T = concat([0], [substvec(q, ['x, 'y], [R[s[1]], R[s[2]]])
		    | s <- duads]);
   M = Mat([[t^i | t <- T]~ | i <- [0..15]]);
   h = 'z * h1;
   red = [polredabs(q, 1) | q <- factor(h)[,1]];
   F = [r[1] | r <- red];
   C = matconcat([algebra_homomorphism_matrix(h, r[1], r[2]) | r <- red]~);
   D = M * C^-1;
   Phi = bestappr(D~ * Winv * D);
   [F, F, Phi];
}

\\ partitions of 5 and 6 leading to a_p = 1 mod 2
ap_1 = [[1, 1, 1, 3], [1, 1, 3], [1, 2, 3], [1, 5], [2, 3], [5]];

ap(h, p) =
{
   if(vecsearch(ap_1, apply(poldegree, factormod(h, p)[,1]~)), 1, 0);
}

compare_traces(h, t) =
{
   my(d = poldisc(h), p);
   for(i = 2, 25,
      p = prime(i);
      if(d % p != 0 && t[i] != ap(h, p),
	 return(0)));
   return(1);
}

dual_pair(h, t) =
{
   my(H);
   if(compare_traces(h, t), dual_pair_GSp4_F2(h),
      H = twin(h); compare_traces(H, t), dual_pair_GSp4_F2(H),
      error("no match for h = ", h, ", t = ", t));
}

match() =
{
   my(reps = apply(x -> strsplit(x, "	"), readstr("GSp4_F2-reps.tsv")), D);
   foreach(reps, r,
      D = dual_pair(Polrev(eval(r[3])), eval(r[2]));
      print(concat([r[1], "	", Str(to_lmfdb_format(D))])));
}
