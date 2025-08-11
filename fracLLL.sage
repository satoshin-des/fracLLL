import random
import time
import matplotlib.pyplot as plt

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
    # B, mu = gso(basis, n)
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

def hf(basis, n):
    factor = basis[0].norm() / abs(det(basis)) ^ (1 / n)
    return factor ^ (1 / n)

def rand_mat(n):
    basis = matrix(ZZ, n, n)
    while True:
        for i in xsrange(n):
            for j in xsrange(n):
                s = random.randint(0, 10)
                if (s & 1) == 0:
                    basis[i, j] = random.randint(100, 10000)
                else:
                    basis[i, j] = -random.randint(100, 10000)
        if basis.det() != 0:
            return basis

if __name__ == '__main__':
    X = []
    frac_lll_time = []
    frac_lll_gso_time = []
    lll_time = []
    for N in xsrange(2, 101, 2):
        X.append(N)
        b = rand_mat(N)
        c = copy(b)
        d = copy(b)
        
        s = time.perf_counter()
        b = frac_lll(b, N, 99, 100)
        e = time.perf_counter()
        frac_lll_gso_time.append(e - s)
        print(f"{N}: fracLLL ended {e - s}")
        s = time.perf_counter()
        c = c.LLL(0.99)
        e = time.perf_counter()
        lll_time.append(e - s)
        print(f"{N}: LLL ended     {e - s}")
        print("-----------------")
    
    _, ax = plt.subplots()
    ax.plot(X, frac_lll_gso_time, marker="", label="fracLLL")
    ax.plot(X, lll_time, marker="", label="LLL")
    
    plt.legend()
    ax.set_xlabel("dimension")
    ax.set_ylabel("run-time[s]")
    plt.savefig(f"time_conp.png")
    plt.show()
