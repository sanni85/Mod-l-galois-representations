/*
  Dual pairs of algebras for 2-dimensional irreducible representations
  of Gal(Qbar/Q) over F_2
*/

W = 1/4 * [1, 1, 1, 1; 1, 1, -1, -1; 1, -1, 1, -1; 1, -1, -1, 1];

\\ f irreducible cubic polynomial (C_3 or S_3)
dual_pair_GL2_F2(f) =
{
   my(g = nfsplitting(subst(f, 'x, 'y)),
      C = matconcat([1, 0; 0, Mat([powers(r, 2)~ | r <- nfroots(g, f)])]),
      F = ['x, f],
      Phi = C * W * C~);
   [F, F, Phi];
}
