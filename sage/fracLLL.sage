NUM = 0
DEN = 1

load("fracGSO.sage")
load("fracSize.sage")

def update_gso_swap(mu: sage.matrix.matrix_integer_dense.Matrix_integer_dense, B: sage.modules.vector_integer_dense.Vector_integer_dense, k: int, n: int) -> tuple:
    """Update Gram-Schmidt orthogonalization informations by exchanging two basis vector b_{k-1} and b_k

    Args:
        mu (sage.matrix.matrix_integer_dense.Matrix_integer_dense): Gram-Schmidt orthogonalization coefficients matrix
        B (sage.modules.vector_integer_dense.Vector_integer_dense): Squared norms of Gram-Schmidt orthogonalization vectors
        k (int): index
        n (int): lattice rank

    Returns:
        tuple: Tuple of Gram-Schmdit orthogonalizartion informations
    """
    C_num = copy(B[NUM])
    C_den = copy(B[DEN])
    nu_num = copy(mu[NUM])
    nu_den = copy(mu[DEN])

    C_num[k - 1] = B[NUM][k] * B[DEN][k - 1] * mu[DEN][k, k - 1] ^ 2 + mu[NUM][k, k - 1] ^ 2 * B[DEN][k] * B[NUM][k - 1]
    C_den[k - 1] = B[DEN][k] * B[DEN][k - 1] * mu[DEN][k, k - 1] ^ 2
    d = gcd(C_num[k - 1], C_den[k - 1])
    C_num[k - 1] //= d
    C_den[k - 1] //= d
    C_num[k] = B[NUM][k - 1] * B[NUM][k] * C_den[k - 1]
    C_den[k] = B[DEN][k - 1] * B[DEN][k] * C_num[k - 1]
    d = gcd(C_num[k], C_den[k])
    C_num[k] //= d
    C_den[k] //= d

    nu_num[k, k - 1] = mu[NUM][k, k - 1] * B[NUM][k - 1] * C_den[k - 1]
    nu_den[k, k - 1] = mu[DEN][k, k - 1] * B[DEN][k - 1] * C_num[k - 1]
    d = gcd(nu_num[k, k - 1], nu_den[k, k - 1])
    nu_num[k, k - 1] //= d
    nu_den[k, k - 1] //= d
    for j in xsrange(k - 1):
        nu_num[k - 1, j] = mu[NUM][k, j]
        nu_den[k - 1, j] = mu[DEN][k, j]
        d = gcd(nu_num[k - 1, j], nu_den[k - 1, j])
        nu_num[k - 1, j] //= d
        nu_den[k - 1, j] //= d
        nu_num[k, j] = mu[NUM][k - 1, j]
        nu_den[k, j] = mu[DEN][k - 1, j]
        d = gcd(nu_num[k, j], nu_den[k, j])
        nu_num[k, j] //= d
        nu_den[k, j] //= d
    for i in xsrange(k + 1, n):
        nu_num[i, k] = mu[NUM][i, k - 1] * mu[DEN][k, k - 1] * mu[DEN][i, k] - mu[NUM][k, k - 1] * mu[NUM][i, k] * mu[DEN][i, k - 1]
        nu_den[i, k] = mu[DEN][i, k - 1] * mu[DEN][k, k - 1] * mu[DEN][i, k]
        d = gcd(nu_num[i, k], nu_den[i, k])
        nu_num[i, k] //= d
        nu_den[i, k] //= d
        nu_num[i, k - 1] = mu[NUM][i, k] * nu_den[k, k - 1] * nu_den[i, k] + mu[DEN][i, k] * nu_num[k, k - 1] * nu_num[i, k]
        nu_den[i, k - 1] = mu[DEN][i, k] * nu_den[k, k - 1] * nu_den[i, k]
        d = gcd(nu_num[i, k - 1], nu_den[i, k - 1])
        nu_num[i, k - 1] //= d
        nu_den[i, k - 1] //= d
    
    return [C_num, C_den], [nu_num, nu_den]

def frac_lll(basis: sage.matrix.matrix_integer_dense.Matrix_integer_dense, n: int, a: int, b: int) -> sage.matrix.matrix_integer_dense.Matrix_integer_dense:
    """Apply fracLLL i.e. LLL-reduction algorithms without FPA

    Args:
        basis (sage.matrix.matrix_integer_dense.Matrix_integer_dense): lattice basis
        n (int): lattice rank
        a (int): numerator of reduction parameter
        b (int): denominator of reduction parameter

    Returns:
        sage.matrix.matrix_integer_dense.Matrix_integer_dense: reduced basis
    """
    B, mu = frac_gso(basis, n)
    k = 1
    while k < n:
        for j in xsrange(k - 1, -1, -1):
            basis, mu = frac_size(basis, mu, k, j)
        
        if B[NUM][k] * B[DEN][k - 1] * b * mu[DEN][k, k - 1] ^ 2 >= (a * mu[DEN][k, k - 1] ^ 2 - b * mu[NUM][k, k - 1] ^ 2) * B[NUM][k - 1] * B[DEN][k]:
            k += 1
        else:
            basis[k], basis[k - 1] = copy(basis[k - 1]), copy(basis[k])
            B, mu = update_gso_swap(mu, B, k, n)
            if k - 1 > 1:
                k -= 1
            else:
                k = 1
    return basis
    