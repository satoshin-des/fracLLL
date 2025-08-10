NUM = 0
DEN = 1

def dij(G, i, j):
    if (i == 0) and (j == 0):
        return 1
    else: 
        d = matrix(ZZ, i, i)
        for k in xsrange(i):
            if k < j:
                d[k] = G[k, : i]
            else:
                d[k] = G[k + 1, : i]
        return d.det()

def frac_gso(basis, n):
    G = basis * basis.transpose()
    B_num = vector(ZZ, n)
    B_den = vector(ZZ, n)
    mu_num = identity_matrix(ZZ, n)
    mu_den = matrix(ZZ, n)
    for i in xsrange(n):
        for j in xsrange(n):
            mu_den[i, j] = 1

    for i in xsrange(n):
        sum = 0
        for j in xsrange(i + 1):
            for k in xsrange(i + 1):
                if ((j + k) & 1) == 0:
                    sum += dij(G, i, j) * dij(G, i, k) * G[j, k]
                else:
                    sum -= dij(G, i, j) * dij(G, i, k) * G[j, k]
        B_num[i] = sum
        B_den[i] = dij(G, i, i) ^ 2

        for j in xsrange(i):
            sum = 0
            for k in xsrange(j + 1):
                if ((j + k) & 1) == 0:
                    sum += dij(G, j, k) * G[i, k]
                else:
                    sum -= dij(G, j, k) * G[i, k]
            mu_num[i, j] = dij(G, j, j) * sum
            
            sum = 0
            for k in xsrange(j + 1):
                for l in xsrange(j + 1):
                    if ((k + l) & 1) == 0:
                        sum += dij(G, j, k) * dij(G, j, l) * G[k, l]
                    else:
                        sum -= dij(G, j, k) * dij(G, j, l) * G[k, l]
            mu_den[i, j] = sum
    return [B_num, B_den], [mu_num, mu_den]

def frac_size(basis, mu, i, j):
    if 2 * abs(mu[NUM][i, j]) > abs(mu[DEN][i, j]):
        q = ZZ(round(mu[NUM][i, j] / mu[DEN][i, j]))
        basis[i] -= q * basis[j]
        for l in xsrange(j + 1):
            mu[NUM][i, l] = mu[NUM][i, l] * mu[DEN][j, l] - q * mu[DEN][i, l] * mu[NUM][j, l]
            mu[DEN][i, l] = mu[DEN][i, l] * mu[DEN][j, l]
            d = gcd(mu[NUM][i, l], mu[DEN][i, l])
            mu[NUM][i, l] //= d
            mu[DEN][i, l] //= d
    return basis, mu

def update_gso_swap(mu, B, k, n):
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

def frac_lll(basis, n, a, b):
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

if __name__ == '__main__':
    N = 20
    b = random_matrix(ZZ, N, N)
    c = copy(b)
    print(frac_lll(b, N, 99, 100))
    print("-----------------")
    print(c.LLL(0.99))
