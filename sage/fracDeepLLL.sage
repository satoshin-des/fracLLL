load("fracGSO.sage")
load("fracSize.sage")

def deep_insertion(basis: sage.matrix.matrix_integer_dense.Matrix_integer_dense, i: int, k: int):
    t = basis[k, :]
    for j in xsrange(k, i, -1):
        basis[j] = basis[j - 1, :]
    basis[i] = t[:]
    return basis

def frac_deep_lll(basis: sage.matrix.matrix_integer_dense.Matrix_integer_dense, n: int, a: int, b: int):
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
                B, mu = frac_gso(basis, n)
                if i > 1:
                    k = i - 1
                else:
                    k = 0
        k += 1
