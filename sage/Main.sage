import random
import time
import matplotlib.pyplot as plt

load("fracLLL.sage")

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
