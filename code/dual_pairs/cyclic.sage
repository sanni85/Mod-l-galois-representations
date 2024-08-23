from dual_pairs.dual_pair_from_cyclic_field import dual_pair_from_cyclic_field

def test(f):
    D = dual_pair_from_cyclic_field(f, GF(2))
    A = D.algebra1()
    B = D.algebra2()
    Phi = D.phi()
    assert all(M == 1 for M in A._basis_matrices())
    assert all(M == 1 for M in B._basis_matrices())
    den = Phi.denominator()
    print([[list(f) for f in A._polys],
           [list(g) for g in B._polys],
           [den, [list(r) for r in den * Phi]]])

R.<x> = QQ[]
test(x^3 - x^2 - 2*x + 1)   # conductor 7^2
test(x^3 - 3*x - 1)         # conductor 3^4
test(x^3 - x^2 - 4*x - 1)   # conductor 13^2
test(x^3 - x^2 - 6*x + 7)   # conductor 19^2
test(x^3 - x^2 - 10*x + 8)  # conductor 31^2
