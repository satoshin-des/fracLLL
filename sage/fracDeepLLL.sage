NUM = 0
DEN = 1

load("fracGSO.sage")
load("fracSize.sage")
load("fracLLL.sage")

def deep_insertion(basis: sage.matrix.matrix_integer_dense.Matrix_integer_dense, i: int, k: int) -> sage.matrix.matrix_integer_dense.Matrix_integer_dense:
    """Apply deep-insertion(i, k) to lattice basis

    Args:
        basis (sage.matrix.matrix_integer_dense.Matrix_integer_dense): lattice basis
        i (int): index
        k (int): index

    Returns:
        sage.matrix.matrix_integer_dense.Matrix_integer_dense: Deep-inserted basis
    """
    t = basis[k, :]
    for j in xsrange(k, i, -1):
        basis[j] = basis[j - 1, :]
    basis[i] = t[:]
    return basis

"""
def update_gso_deep_insertion(mu: list, B: list, G: ZZ, i: int, k: int, n: int) -> tuple:
    C_num = copy(B[NUM])
    C_den = copy(B[DEN])
    nu_num = copy(mu[NUM])
    nu_den = copy(mu[DEN])
    D_num = vector(ZZ, k)
    D_den = vector(ZZ, k)
    P = vector(ZZ, n)

    P[0] = 1
    for j in xsrange(1, n):
        P[j] = P[j - 1] * mu[DEN][k, j - 1] ^ 2 * B[DEN][j - 1]
    
    for l in xsrange(n):
        D_den[l] = P[l]
        sum = 0
        for j in xsrange(l):
            sum += mu[NUM][k, j] * B[NUM][j] * prod // (mu[DEN][k, j] ^ 2 * B[DEN][j])
        D_num[l] = G * P[l] - sum
    
    C_num[i] = D_num[i]
    C_den[i] = D_den[i]
    for j in xsrange(i + 1, k + 1):
        C_num[j] = D_den[j - 1] * D_num[j] * B[NUM][j - 1]
        C_den[j] = D_num[j - 1] * D_den[j] * B[DEN][j - 1]

    for j in xsrange(i + 1, k + 1):
        pass
"""

def update_gso_conti_swap(mu: list, B: list, i: int, k: int, n: int) -> tuple:
    C_num = copy(B[NUM])
    C_den = copy(B[DEN])
    nu_num = copy(mu[NUM])
    nu_den = copy(mu[DEN])

    for l in xsrange(k, i, -1):
        C_num[l - 1] = B[NUM][l] * B[DEN][l - 1] * mu[DEN][l, l - 1] ^ 2 + mu[NUM][l, l - 1] ^ 2 * B[DEN][l] * B[NUM][l - 1]
        C_den[l - 1] = B[DEN][l] * B[DEN][l - 1] * mu[DEN][l, l - 1] ^ 2
        d = gcd(C_num[l - 1], C_den[l - 1])
        C_num[l - 1] //= d
        C_den[l - 1] //= d
        C_num[l] = B[NUM][l - 1] * B[NUM][l] * C_den[l - 1]
        C_den[l] = B[DEN][l - 1] * B[DEN][l] * C_num[l - 1]
        d = gcd(C_num[l], C_den[l])
        C_num[l] //= d
        C_den[l] //= d

        nu_num[l, l - 1] = mu[NUM][l, l - 1] * B[NUM][l - 1] * C_den[l - 1]
        nu_den[l, l - 1] = mu[DEN][l, l - 1] * B[DEN][l - 1] * C_num[l - 1]
        d = gcd(nu_num[l, l - 1], nu_den[l, l - 1])
        nu_num[l, l - 1] //= d
        nu_den[l, l - 1] //= d
        for j in xsrange(l - 1):
            nu_num[l - 1, j] = mu[NUM][l, j]
            nu_den[l - 1, j] = mu[DEN][l, j]
            d = gcd(nu_num[l - 1, j], nu_den[l - 1, j])
            nu_num[l - 1, j] //= d
            nu_den[l - 1, j] //= d
            nu_num[l, j] = mu[NUM][l - 1, j]
            nu_den[l, j] = mu[DEN][l - 1, j]
            d = gcd(nu_num[l, j], nu_den[l, j])
            nu_num[l, j] //= d
            nu_den[l, j] //= d
        for i in xsrange(l + 1, n):
            nu_num[i, l] = mu[NUM][i, l - 1] * mu[DEN][l, l - 1] * mu[DEN][i, l] - mu[NUM][l, l - 1] * mu[NUM][i, l] * mu[DEN][i, l - 1]
            nu_den[i, l] = mu[DEN][i, l - 1] * mu[DEN][l, l - 1] * mu[DEN][i, l]
            d = gcd(nu_num[i, l], nu_den[i, l])
            nu_num[i, l] //= d
            nu_den[i, l] //= d
            nu_num[i, l - 1] = mu[NUM][i, l] * nu_den[l, l - 1] * nu_den[i, l] + mu[DEN][i, l] * nu_num[l, l - 1] * nu_num[i, l]
            nu_den[i, l - 1] = mu[DEN][i, l] * nu_den[l, l - 1] * nu_den[i, l]
            d = gcd(nu_num[i, l - 1], nu_den[i, l - 1])
            nu_num[i, l - 1] //= d
            nu_den[i, l - 1] //= d
        
        B[NUM] = C_num[:]
        B[DEN] = C_den[:]
        mu[NUM] = nu_num[:, :]
        mu[DEN] = nu_den[:, :]
    
    return [C_num, C_den], [nu_num, nu_den]

def frac_deep_lll(basis: sage.matrix.matrix_integer_dense.Matrix_integer_dense, n: int, a: int, b: int) -> sage.matrix.matrix_integer_dense.Matrix_integer_dense:
    """Compute DeepLLL basis without FPA

    Args:
        basis (sage.matrix.matrix_integer_dense.Matrix_integer_dense): lattice basis
        n (int): lattice rank
        a (int): numerator of reduction parameter
        b (int): denominator of reduction parameter

    Returns:
        sage.matrix.matrix_integer_dense.Matrix_integer_dense: DeepLLL-reduced basis
    """
    B, mu = frac_gso(basis, n)

    k = 1
    while k < n:
        for j in xsrange(k - 1, -1, -1):
            basis, mu = frac_size(basis, mu, k, j)
        
        i = 0
        prod = 1
        sum = 0
        g = basis[k].inner_product(basis[k])
        while i < k:
            if (b * g * B[DEN][i] - a * B[NUM][i]) * prod >= b * B[DEN][i] * sum:
                prod *= mu[DEN][k, i] ^ 2 * B[DEN][i]
                i += 1
                sum = 0
                for j in xsrange(i):
                    sum += mu[NUM][k, j] * B[NUM][j] * prod // (mu[DEN][k, j] ^ 2 * B[DEN][j])
            else:
                basis = deep_insertion(basis, i, k)
                B, mu = update_gso_conti_swap(mu, B, i, k, n)
                # for j in xsrange(k, i, -1):
                #     B, mu = update_gso_swap(mu, B, j, n)
                
                if i > 1:
                    k = i - 1
                else:
                    k = 0
        k += 1
    return basis
