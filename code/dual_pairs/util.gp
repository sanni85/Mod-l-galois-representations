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

to_lmfdb_format(D) =
{
   my(A, B, Phi, d);
   [A, B, Phi] = D;
   d = denominator(Phi);
   [apply(Vecrev, A), apply(Vecrev, B), [d, [v~ | v <- d*Phi~]]];
}
