NUM = 0
DEN = 1

def frac_size(basis: sage.matrix.matrix_integer_dense.Matrix_integer_dense, mu: list, i: int, j: int) -> tuple:
    """Size-reduction with rational number representation.

    Args:
        basis (sage.matrix.matrix_integer_dense.Matrix_integer_dense): Lattice basis.
        mu (list): Gram-Schmidt orthogonalization coefficients matrix
        i (int): index
        j (int): index

    Returns:
        tuple: Tuple of updated lattice basis and  Gram-Schmidt orthogonalization coefficients matrix
    """
    if 2 * abs(mu[NUM][i, j]) > abs(mu[DEN][i, j]):
        q = ZZ(round(mu[NUM][i, j] / mu[DEN][i, j]))
        basis[i] -= q * basis[j]
        for l in xsrange(j + 1):
            mu[NUM][i, l] = mu[NUM][i, l] * mu[DEN][j, l] - q * mu[DEN][i, l] * mu[NUM][j, l]
            mu[DEN][i, l] *= mu[DEN][j, l]
            d = gcd(mu[NUM][i, l], mu[DEN][i, l])
            mu[NUM][i, l] //= d
            mu[DEN][i, l] //= d
    return basis, mu
