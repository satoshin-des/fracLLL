import numpy as np

load("fracDeepLLL.sage")

def update_gso(i, k, B, mu, n):
    P = np.zeros(n)
    D = np.zeros(n)
    S = np.zeros(n)

    P[k] = D[k] = B[k]
    for j in xsrange(k - 1, i - 1, -1):
        P[j] = mu[k, j] * B[j]
        D[j] = D[j + 1] + mu[k, j] * P[j]
    
    for j in xsrange(k, i, -1):
        t = mu[k, j - 1] / D[j]
        for l in xsrange(n - 1, k, -1):
            S[l] += mu[l, j] * P[j]
            mu[l, j] = mu[l, j - 1] - t * S[l]
        for l in xsrange(k, j, -1):
            S[l] += mu[l - 1, j] * P[j]
            mu[l, j] = mu[l - 1, j - 1] - t * S[l]
    
    t = 1. / D[i]

    for l in xsrange(n - 1, k, -1):
        mu[l, i] = t * (S[l] + mu[l, i] * P[i])
    for l in xsrange(k, i + 1, -1):
        mu[l, i] = t * (S[l] + mu[l - 1, i] * P[i])
    
    mu[i + 1, i] = t * P[i]
    for j in xsrange(i):
        eps = mu[k, j]
        for l in xsrange(k, i, -1):
            mu[l, j] = mu[l - 1, j]
        mu = eps
    
    for j in xsrange(k, i, -1):
        B[j] = D[j] * B[j - 1] / D[j - 1]
    B[i] = D[i]

    return B, mu

def deep_lll(basis, delta, n):
    b, mu = basis.gram_schmidt()

    mu = np.zeros((n, n))
    B = np.zeros(n)
    for i in xsrange(n):
        B[i] = b.norm() ^ 2

    k = 1
    while k < n:
        for j in xsrange(k - 1, -1, -1):
            if mu[k, j] > 0.5 or mu[k, j] < -0.5:
                q = round(mu[k, j])
                basis[k] -= q * basis[j]
                mu[k, : j + 1] -= q * mu[j, : j + 1]
            
        C = basis[k].norm() ^ 2
        i = 0
        while i < k:
            if i < k:
                C -= mu[k, i] ^ 2 * B[i]
                i += 1
            else:
                basis = deep_insertion(basis, i, k)
                B, mu = update_gso(i, k, B, mu, n)
                if i > 1:
                    k = i - 1
                else:
                    k = 0
        k += 1

    return basis
