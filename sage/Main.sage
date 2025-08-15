import random
import time
import matplotlib.pyplot as plt

load("fracLLL.sage")
load("fracDeepLLL.sage")
load("DeepLLL.sage")

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

def rhf(basis: sage.matrix.matrix_integer_dense.Matrix_integer_dense, n: int) -> RR:
    hermite_factor = basis[0].norm() / abs(basis.det()) ^ (1 / n)
    return RR(hermite_factor ^ (1 / n))

def sl(basis: sage.matrix.matrix_integer_dense.Matrix_integer_dense, n: int) -> RR:
    S = 0.
    T = 0.
    B, _ = basis.gram_schmidt()
    for i in xsrange(n):
        S += RR((i + 1) * 2 * log(B[i].norm()))
        T += 2 * log(B[i].norm())
    return RR(abs(12 * (S - 0.5 * (n + 1) * T) / (n * (n * n - 1))))

def collect_time(red: function) -> None:
    X = []
    frac_red_time = []
    lll_time = []
    for N in xsrange(2, 101, 2):
        X.append(N)
        b = rand_mat(N)
        c = copy(b)
        d = copy(b)
        
        s = time.perf_counter()
        b = red(b, N, 99, 100)
        e = time.perf_counter()
        frac_red_time.append(e - s)
        print(f"{N}: fracLLL ended {e - s}")
        s = time.perf_counter()
        c = c.LLL(0.99)
        e = time.perf_counter()
        lll_time.append(e - s)
        print(f"{N}: LLL ended     {e - s}")
        print("-----------------")
    
    _, ax = plt.subplots()
    ax.plot(X, frac_red_time, marker="", label="fracLLL")
    ax.plot(X, lll_time, marker="", label="LLL")
    
    plt.legend()
    ax.set_xlabel("dimension")
    ax.set_ylabel("run-time[s]")
    plt.savefig(f"time_conp.png")
    plt.show()

def collect_rhf(red: function) -> None:
    X = []
    frac_red_rhf = []
    lll_rhf = []
    for N in xsrange(40, 151, 10):
        X.append(N)
        b = deep_lll(rand_mat(N).LLL(0.99), 0.99, N)
        c = copy(b)
        d = copy(b)
        
        b = red(c, N, 99, 100)
        frac_red_rhf.append(rhf(b, N))
        print(f"{N}: frac ended {frac_red_rhf[-1]}")

        # b = deep_lll(d, 0.99, N)
        lll_rhf.append(rhf(d, N))
        print(f"{N}: LLL ended  {lll_rhf[-1]}")
        print("-----------------")
    
    _, ax = plt.subplots()
    ax.plot(X, frac_red_rhf, marker="", label="fracLLL")
    ax.plot(X, lll_rhf, marker="", label="LLL")
    
    plt.legend()
    ax.set_xlabel("dimension")
    ax.set_ylabel("RHF")
    plt.savefig(f"rhf_conp.png")
    plt.show()

def collect_sl(red: function) -> None:
    X = []
    frac_red_sl = []
    red_sl = []
    for N in xsrange(40, 151, 10):
        X.append(N)
        b = deep_lll(rand_mat(N).LLL(0.99), 0.99, N)
        c = copy(b)
        d = copy(b)
        
        b = red(c, N, 99, 100)
        frac_red_sl.append(sl(b, N))
        print(f"{N}: frac ended {frac_red_sl[-1]}")

        # b = deep_lll(d, 0.99, N)
        red_sl.append(sl(d, N))
        print(f"{N}: LLL ended  {red_sl[-1]}")
        print("-----------------")
    
    _, ax = plt.subplots()
    ax.plot(X, frac_red_sl, marker="", label="fracLLL")
    ax.plot(X, red_sl, marker="", label="LLL")
    
    plt.legend()
    ax.set_xlabel("dimension")
    ax.set_ylabel("GSA-slope")
    plt.savefig(f"sl_conp.png")
    plt.show()

def collect_rhf_and_sl(red: function):
    X = []
    frac_red_sl = []
    red_sl = []
    frac_red_rhf = []
    lll_rhf = []
    for N in xsrange(40, 151, 10):
        X.append(N)
        b = deep_lll(rand_mat(N).LLL(0.99), 0.99, N)
        c = copy(b)
        d = copy(b)
        
        b = red(c, N, 99, 100)
        frac_red_sl.append(sl(b, N))
        print(f"{N}: frac ended {frac_red_sl[-1]}")
        frac_red_rhf.append(rhf(b, N))
        print(f"{N}: frac ended {frac_red_rhf[-1]}")

        # b = deep_lll(d, 0.99, N)
        red_sl.append(sl(d, N))
        print(f"{N}: LLL ended  {red_sl[-1]}")
        lll_rhf.append(rhf(d, N))
        print(f"{N}: LLL ended  {lll_rhf[-1]}")
        print("-----------------")
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(X, frac_red_sl, marker="", label="fracLLL", color='b')
    ax1.plot(X, red_sl, marker="", label="LLL", color='c')
    ax1.set_xlabel("dimension")
    ax1.set_ylabel("GSA-slope")
    
    ax2 = ax1.twinx()
    ax2.plot(X, frac_red_rhf, marker="", label="fracLLL")
    ax2.plot(X, lll_rhf, marker="", label="LLL")
    ax2.set_ylabel("RHF")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2)

    plt.tick_params()

    plt.savefig(f"sl_rhf_conp.png")
    plt.show()

if __name__ == '__main__':
    #collect_rhf(frac_deep_lll)
    #collect_sl(frac_deep_lll)
    # collect_time(frac_deep_lll)
    collect_rhf_and_sl(frac_deep_lll)
