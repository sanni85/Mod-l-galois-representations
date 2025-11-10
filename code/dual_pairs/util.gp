/*
  Some utility functions
*/

find(v, P) = for(j = 1, #v, if(P(v[j]), return(j)));

\\ vector p such that vecextract(x, p) = y
findvec(x, y) = [find(x, t -> t == z) | z <- y];

\\ block matrix formed of square matrices
matblock(M) =
{
   matconcat(matdiagonal([m | m <- M, m != [;]]));
}

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

\\ matrix of the algebra endomorphism Q[x]/(f) -> Q[y]/(g)
\\ defined by sending x to a (assuming g | f(a))
algebra_homomorphism_matrix(f, g, a) =
{
   my(d = poldegree(f), e = poldegree(g));
   Mat([Colrev(liftpol(b), e) | b <- powers(Mod(a, g), d - 1)]);
}

nf_zk_matrix(K) =
{
   my(d = poldegree(K.pol));
   Mat(apply(x -> Colrev(x, d), K.zk));
}

nf_homomorphism_matrix(K, L, a) =
{
   L[8] * algebra_homomorphism_matrix(K.pol, L.pol, a) * nf_zk_matrix(K);
}

\\ express an element of L = K[x]/(g) on the basis
\\ (w_1, ..., w_n, x, w_1*x, ..., w_n*x^(deg(g) - 1))
\\ where (w_1, ..., w_n) is an integral basis of K
algtobasis_rel(K, g, a) =
{
   a = liftpol(Mod(a, g));
   concat([nfalgtobasis(K, polcoef(a, i, 'x)) | i <- [0..poldegree(g) - 1]]);
}

homomorphism_matrix_rel(K, g, a) =
{
   Mat([algtobasis_rel(K, g, b) |
	b <- powers(a, poldegree(K.pol) - 1)]) * nf_zk_matrix(K);
}

to_lmfdb_format(D) =
{
   my([A, B, Phi] = D,
      d = denominator(Phi));
   [apply(Vecrev, A), apply(Vecrev, B), [d, [v~ | v <- d*Phi~]]];
}

from_lmfdb_format(D) =
{
   my([A, B, Phi] = D);
   [apply(Polrev, A), apply(Polrev, B), 1/Phi[1] * Mat([v~ | v <- Phi[2]])~];
}
