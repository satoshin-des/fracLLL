NUM = 0
DEN = 1

def frac_gso(basis: sage.matrix.matrix_integer_dense.Matrix_integer_dense, n: int) -> tuple:
    """Compute Gram-Schmidt orthogonalization information as rational number

    Args:
        basis (sage.matrix.matrix_integer_dense.Matrix_integer_dense): Lattice basis
        n (int): Rank of lattice

    Returns:
        tuple: Tuple of Gram-Schmidt orthogonalization informations
    """
    G = basis * basis.transpose()
    B_num = vector(ZZ, n)
    B_den = vector(ZZ, n)
    mu_num = identity_matrix(ZZ, n)
    mu_den = matrix(ZZ, n, n)
    D = matrix(ZZ, n, n)
    for i in xsrange(n):
        for j in xsrange(n):
            mu_den[i, j] = 1

    for i in xsrange(n):
        for j in xsrange(n):
            if (i == 0) and (j == 0):
                D[i, j] = 1
            else:
                d = matrix(ZZ, i, i)
                for k in xsrange(i):
                    if k < j:
                        d[k] = G[k, : i]
                    else:
                        d[k] = G[k + 1, : i]
                D[i, j] = d.det()

    for i in xsrange(n):
        sum = 0
        for j in xsrange(i + 1):
            for k in xsrange(i + 1):
                if ((j + k) & 1) == 0:
                    sum += D[i, j] * D[i, k] * G[j, k]
                else:
                    sum -= D[i, j] * D[i, k] * G[j, k]
        B_num[i] = sum
        B_den[i] = D[i, i] * D[i, i]
        d = gcd(B_num[i], B_den[i])
        B_num[i] //= d
        B_den[i] //= d

        for j in xsrange(i):
            sum = 0
            for k in xsrange(j + 1):
                if ((j + k) & 1) == 0:
                    sum += D[j, k] * G[i, k]
                else:
                    sum -= D[j, k] * G[i, k]
            mu_num[i, j] = D[j, j] * sum
            
            sum = 0
            for k in xsrange(j + 1):
                for l in xsrange(j + 1):
                    if ((k + l) & 1) == 0:
                        sum += D[j, k] * D[j, l] * G[k, l]
                    else:
                        sum -= D[j, k] * D[j, l] * G[k, l]
            mu_den[i, j] = sum
            d = gcd(mu_num[i, j], mu_den[i, j])
            mu_num[i, j] //= d
            mu_den[i, j] //= d
    return [B_num, B_den], [mu_num, mu_den]
