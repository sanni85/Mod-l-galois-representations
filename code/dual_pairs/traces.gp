/*
  Frobenius traces for dual pairs of algebras
*/

\r util.gp

dual_pair_import(f) =
{
   my(D = read(f),
      alg1 = [nfinit(Polrev(v)) | v <- D[1]],
      alg2 = [nfinit(Polrev(v)) | v <- D[2]],
      B1 = matblock([K[8] | K <- alg1]),
      B2 = matblock([K[8] | K <- alg2]),
      phi = B1~^-1 * (1/D[3][1] * Mat([v~ | v <- D[3][2]])~) * B2^-1);
   [alg1, alg2, phi];
}

to_col(x) =
{
   if(type(x) == "t_INTMOD", Col(x),
      Mod(Colrev(x.pol, poldegree(x.mod)), x.p));
}

nf_reduce(K, p) =
{
   my(modpr = [nfmodprinit(K, q) | q <- idealprimedec(K, p)]);
   [[nfmodpr(K, 0, q) | q <- modpr],
    Mat([concat([to_col(nfmodpr(K, x, q)) | q <- modpr]) | x <- K.zk])];
}

algebra_reduce(A, p) =
{
   my(T = [nf_reduce(K, p) | K <- A]);
   [concat([t[1] | t <- T]),
    Mod(matblock([t[2] | t <- T]), p)];
}

dual_pair_reduce(D, p) =
{
   my([A1, B1] = algebra_reduce(D[1], p),
      [A2, B2] = algebra_reduce(D[2], p));
   [A1, A2, B1~^-1 * D[3] * B2^-1];
}

\\ A: product of finite fields
splitting_field_degree(A) =
{
   lcm([poldegree(x.mod) | x <- A, type(x) == "t_FFELT"]);
}

splitting_field(D) =
{
   my(p = characteristic(D[3]));
   ffgen([p, lcm(splitting_field_degree(D[1]),
		 splitting_field_degree(D[2]))]);
}

\\ assume deg F0 | deg F1 and F1 a t_FFELT
ff_morphisms(F0, F1) =
{
   if(type(F0) == "t_INTMOD", return(Col(F1^0)));
   my(p = F0.p,
      d = poldegree(F0.mod),
      emb = vectorv(d),
      f = ffembed(F0, F1));
   emb[1] = ffmap(f, powers(ffgen(F0), d - 1));
   for(i = 1, d - 1,
      emb[i + 1] = apply(z -> z^p, emb[i]));
   matconcat(emb);
}

algebra_ff_morphisms(A, F) =
{
   matblock([ff_morphisms(G, F) | G <- A]);
}

group_structure(A1, A2, th, F) =
{
   my(emb1, emb2, n, q, o, z, T, B, C,
      u, v, d, U, elements, j, E, V);

   emb1 = algebra_ff_morphisms(A1, F);
   emb2 = algebra_ff_morphisms(A2, F);
   n = matsize(emb1)[1];
   q = F.p ^ poldegree(F.mod);
   o = gcd(q - 1, n);
   z = polrootsmod(polcyclo(o), F)[1];
   T = apply(w -> fflog(w, z, o), emb1 * th * emb2~) / o;

   \\ Compute integral vectors in the image of T...
   B = matrixqz(T, -1);
   \\ ... and integral vectors mapping to these.
   C = matrixqz(concat(matinverseimage(T, B), matker(T)), -1);

   [u, v, d] = matsnf(C, 5);
   U = apply(frac, T*C*v*d^-1);

   d = vector(#d, i, d[i, i]);
   elements = vector(n);
   j = 0;
   forvec(i = [[0, o - 1] | o <- d],
	  elements[j++] = i~);

   \\ matrix of the pairing with respect to the given elements
   E = matrix(n, n, i, j,
	      frac(sum(k = 1, #d, elements[i][k] * elements[j][k] / d[k])));

   V = apply(frac, Mat([U * e | e <- elements]));
   [d, [vectorv(#d, i, Mod(e[i], d[i])) | e <- elements],
    vecextract(Col(emb1), findvec(Col(V), Col(E)))];
}

frobenius_matrix(D) =
{
   my([A1, A2, Phi] = D,
      p = characteristic(D[3]),
      F = splitting_field(D),
      [d, V, P] = group_structure(A1, A2, Phi~^-1, F),
      basis = vecextract(P, findvec(V, Vec(matid(#d)))),
      basis1 = [[x^p | x <- v] | v <- basis]);
   Mat(vecextract(V, findvec(P, basis1)));
}

frobenius_traces(D) =
{
   my(n = matsize(D[3])[1],
      disc = [K.disc | K <- D[1]]);
   [if(n % p == 0 || find(disc, d -> d % p == 0), -1,
       liftint(trace(frobenius_matrix(dual_pair_reduce(D, p)))))
    | p <- primes(25)];
}
