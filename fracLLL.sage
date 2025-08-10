def dij(G, i, j):
    if (i == 0) and (j == 0):
        return 1
    else: 
        d = matrix(ZZ, i, i)
        for k in xsrange(i):
            if k < j:
                d[k] = G[k, : i]
            elif k > j:
                d[k] = G[k + 1, : i]
        return d.det()

def fracGSO(b, n):
    G = b * b.transpose()
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
            for k in xsrange(i + 1):
                for l in xsrange(j + 1):
                    if ((k + l) & 1) == 0:
                        sum += dij(G, j, k) * dij(G, j, l) * G[k, l]
                    else:
                        sum -= dij(G, j, k) * dij(G, j, l) * G[k, l]
            mu_den[i, j] = sum
    return [B_num, B_den], [mu_num, mu_den]

if __name__ == '__main__':
    N = 10
    b = random_matrix(ZZ, N, N)
    B = vector(QQ, N)
    mu = matrix(QQ, N, N)
    frac_B, frac_mu = fracGSO(b, N)

    for i in xsrange(N):
        for j in xsrange(N):
            mu[i, j] = frac_mu[0][i, j] / frac_mu[1][i, j]
        B[i] = frac_B[0][i] / frac_B[1][i]
    print(B)
    print(mu)
    
    B, mu = b.gram_schmidt()
    # print(B)
    print(mu)
