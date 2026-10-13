import argparse
import itertools
import sys
import numpy as np
from math import sqrt, pi, exp
from scipy.special import erfc
import harm3 as H
import harmn as N


def nearset(d, L, K=1):
    out = []
    for R in itertools.product(range(-K - 3, K + 4), repeat=3):
        y = np.array(d) + np.array(R) * np.asarray(L)
        if max(abs(y)) <= K:
            out.append(y)
    return out


def latfar(d, L, J, a=None, tol=1e-17, K=1):
    L = np.asarray(L, float)
    V = L.prod()
    if a is None:
        a = 1.0
    x = -np.asarray(d, float)
    near = nearset(d, L, K)
    nearkey = {tuple(int(v) for v in y) for y in near}
    rc = (sqrt(-np.log(tol)) + sqrt(J) + 3) / a
    nmax = [int(np.ceil(rc / L[k])) + 2 for k in range(3)]
    n, m = H.nm(J)
    Jn = n
    RR = np.array(list(itertools.product(*[range(-nmax[k], nmax[k] + 1) for k in range(3)])), float)
    off = np.rint(np.array(d)[None, :] + RR * L[None, :]).astype(int)
    isnear = np.array([tuple(o) in nearkey for o in off.tolist()])
    Y = x[None, :] - RR * L[None, :]
    r = np.sqrt((Y * Y).sum(1))
    sel = (~isnear) & (r <= rc + 3)
    Y, r = Y[sel], r[sel]
    Bt = np.zeros((J + 1, len(r)))
    Bt[0] = erfc(a * r) / r
    g = 1 / (a * sqrt(pi))
    e = np.exp(-a * a * r * r)
    for j in range(1, J + 1):
        g = g * 2 * a * a / (2 * j - 1)
        Bt[j] = (Bt[j - 1] + g * e) / (r * r)
    S = np.zeros(H.nterm(J), complex)
    for i in range(0, len(r), 4096):
        S += (N.reg(Y[i:i + 4096], J) * Bt[Jn][:, i:i + 4096].T).sum(0)
    kc = 2 * a * (sqrt(-np.log(tol)) + sqrt(J) + 3)
    kmax = [int(np.ceil(kc * L[k] / (2 * pi))) + 2 for k in range(3)]
    K = np.array(list(itertools.product(*[np.arange(-kmax[k], kmax[k] + 1) for k in range(3)])), float)
    K = K[(K != 0).any(1)] * 2 * pi / L
    k2 = (K * K).sum(1)
    w = 4 * pi / V * np.exp(-k2 / (4 * a * a)) / k2 * np.exp(1j * (K @ x))
    keep = np.abs(w) > 1e-300
    K, w, k2 = K[keep], w[keep], k2[keep]
    km = np.sqrt(k2).max()
    sc = np.ones(J + 1)
    for j in range(1, J + 1):
        sc[j] = sc[j - 1] * km / (2 * j - 1)
    rec = np.zeros(H.nterm(J), complex)
    for i in range(0, len(K), 2048):
        rec += (w[i:i + 2048, None] * N.reg(K[i:i + 2048], J, km)).sum(0)
    S += rec * sc[Jn] * (1j)**Jn * (-1.0)**Jn
    for y0 in near:
        y = -np.asarray(y0, float)
        u = (y * y).sum()
        Ej = np.zeros(J + 1)
        pre = 2 * a / sqrt(pi)
        for j in range(J + 1):
            if j > 0:
                pre = pre * 2 * a * a / (2 * j - 1)
            s, tk, i = 0.0, 1.0, 0
            while True:
                term = tk / (2 * j + 2 * i + 1)
                s += term
                i += 1
                tk = tk * (-a * a * u) / i
                if abs(tk) < 1e-18 * max(abs(s), 1e-300) and i > 3:
                    break
            Ej[j] = pre * s
        S -= N.reg(y, J) * Ej[Jn]
    S[0] -= pi / (V * a * a)
    S[0] -= 2 * pi / (3 * V) * (x * x).sum()
    R1 = N.reg(x, 1)
    for mm in (-1, 0, 1):
        S[H.ix(1, mm)] += 4 * pi / (3 * V) * R1[H.ix(1, mm)]
    return S


def farops(cells, p, fe, K=1):
    n1, n2, n3 = cells
    L = np.array([n1, n2, n3], float)
    nn, mm = H.nm(p)
    nh = [(a, b) for a in range(p + 1) for b in range(a + 1)]
    out = {}
    for d in itertools.product(range(-(n1 - 1), n1), range(-(n2 - 1), n2), range(-(n3 - 1), n3)):
        S = latfar(d, L, 2 * p, K=K)
        T = np.zeros((len(nh), len(nn)), complex)
        for i, (n, m) in enumerate(nh):
            for j, (k, l) in enumerate(zip(nn, mm)):
                T[i, j] = (-1)**abs(l) * fe[H.ix(n, m), j] * S[H.ix(n + k, m - l)]
        out[d] = T
    return out


if __name__ == "__main__":
    pa = argparse.ArgumentParser(prog="periodicn")
    pa.add_argument("-c", nargs=3, type=int, default=[1, 1, 1])
    pa.add_argument("-p", type=int, default=10)
    pa.add_argument("-T", required=True)
    pa.add_argument("-o", required=True)
    pa.add_argument("-K", type=int, default=1)
    a = pa.parse_args()
    if a.p < 1 or min(a.c) < 1:
        sys.stderr.write("periodicn: error: need -p >= 1 and positive cell counts\n")
        sys.exit(1)
    with open(a.T, "rb") as f:
        tp = np.fromfile(f, np.int64, 1)[0]
        if tp != a.p:
            sys.stderr.write("periodicn: error: table '%s' is for p = %d\n" % (a.T, tp))
            sys.exit(1)
        for n in range(a.p + 1):
            np.fromfile(f, float, (2 * n + 1)**2)
        np.fromfile(f, float, 3 * (a.p + 1)**3)
        nt = (a.p + 1)**2
        fe = np.fromfile(f, float, nt * nt).reshape(nt, nt)
    ops = farops(a.c, a.p, fe, a.K)
    with open(a.o, "wb") as f:
        np.array([a.c[0], a.c[1], a.c[2], a.p, a.K], np.int64).tofile(f)
        for d in itertools.product(range(-(a.c[0] - 1), a.c[0]), range(-(a.c[1] - 1), a.c[1]), range(-(a.c[2] - 1), a.c[2])):
            ops[d].astype(np.complex128).tofile(f)
