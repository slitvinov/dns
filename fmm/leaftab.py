import argparse
import sys
import numpy as np
import mpmath as mp

pa = argparse.ArgumentParser(prog="leaftab")
pa.add_argument("-B", type=int, required=True)
pa.add_argument("-p", type=int, required=True)
pa.add_argument("-o", required=True)
a = pa.parse_args()
B, p = a.B, a.p
if B < 1 or B & (B - 1) or p < 1 or p > 40:
    sys.stderr.write("leaftab: error: need -B a power of 2, 1 <= -p <= 40\n")
    sys.exit(1)
mp.mp.dps = 40
nt = (p + 1) ** 2
B3 = B ** 3


def ix(n, m):
    return n * n + n + m


def lc(i):
    x = y = z = 0
    b = 0
    while (1 << b) < B:
        x |= ((i >> (3 * b + 2)) & 1) << b
        y |= ((i >> (3 * b + 1)) & 1) << b
        z |= ((i >> (3 * b)) & 1) << b
        b += 1
    return x, y, z


def reg(X, Y, Z):
    r2 = X * X + Y * Y + Z * Z
    w = mp.mpc(X, Y)
    R = {}
    um = mp.mpc(1)
    for m in range(p + 1):
        if m > 0:
            um = um * w * mp.sqrt(mp.mpf(2 * m - 1) / (2 * m))
        R[(m, m)] = um
        u1, u2 = um, mp.mpc(0)
        for n in range(m + 1, p + 1):
            u = ((2 * n - 1) * Z * u1 - mp.sqrt((n - 1) ** 2 - m * m) * r2 * u2) / mp.sqrt(n * n - m * m)
            R[(n, m)] = u
            u2, u1 = u1, u
    F = {}
    for (n, m), v in R.items():
        F[(n, m)] = v
        if m > 0:
            F[(n, -m)] = (-1) ** m * mp.conj(v)
    return F


P = np.zeros((B3, 2 * nt))
V = np.zeros((6 * nt, 3 * B3))
pi4 = 4 * mp.pi
for i in range(B3):
    X = [mp.mpf(2 * c + 1 - B) / (2 * B) for c in lc(i)]
    F = reg(*X)
    for (n, m), v in F.items():
        P[i, 2 * ix(n, m)] = float(v.real)
        P[i, 2 * ix(n, m) + 1] = float(v.imag)
    F = reg(*[-x for x in X])
    W = [dict(), dict(), dict()]

    def put(d, key, v):
        W[d][key] = W[d].get(key, mp.mpc(0)) + v

    for n in range(p):
        for mm in range(-n, n + 1):
            cu = mp.conj(F[(n, mm)])
            ap = mp.sqrt((n + mm + 1) * (n + mm + 2))
            am = mp.sqrt((n - mm + 1) * (n - mm + 2))
            put(0, (n + 1, mm + 1), -cu * ap / 2)
            put(0, (n + 1, mm - 1), cu * am / 2)
            put(1, (n + 1, mm + 1), 1j * cu * ap / 2)
            put(1, (n + 1, mm - 1), 1j * cu * am / 2)
            put(2, (n + 1, mm), -mp.sqrt((n + 1) ** 2 - mm * mm) * cu)
    for e, terms in enumerate((((2, 1, 1), (1, 2, -1)), ((0, 2, 1), (2, 0, -1)), ((1, 0, 1), (0, 1, -1)))):
        for c, d, s in terms:
            for (n, m), w in W[d].items():
                col = c * 2 * nt + 2 * ix(n, m)
                V[col, 3 * i + e] = float(s * w.real / pi4)
                V[col + 1, 3 * i + e] = float(-s * w.imag / pi4)
with open(a.o, "wb") as f:
    np.array([B, p], np.int64).tofile(f)
    P.tofile(f)
    V.tofile(f)
