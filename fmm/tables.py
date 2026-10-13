import argparse
import sys
import numpy as np
import mpmath as mp
from math import comb
import harm3 as H

pa = argparse.ArgumentParser(prog="tables")
pa.add_argument("-p", type=int, required=True)
pa.add_argument("-o", required=True)
a = pa.parse_args()
p = a.p
if p < 1 or p > 40:
    sys.stderr.write("tables: error: need 1 <= -p <= 40\n")
    sys.exit(1)
mp.mp.dps = 50
out = [np.array([p], np.int64).tobytes()]
for n in range(p + 1):
    w = 2 * n + 1
    A = mp.zeros(w, w)
    for i in range(w):
        m = i - n
        if m < n:
            A[i + 1, i] = -mp.pi / 4 * mp.sqrt(n * (n + 1) - m * (m + 1))
        if m > -n:
            A[i - 1, i] = mp.pi / 4 * mp.sqrt(n * (n + 1) - m * (m - 1))
    E = mp.expm(A)
    out.append(np.array([[float(E[i, j]) for j in range(w)] for i in range(w)], float).tobytes())
P1 = p + 1
fa = np.zeros((P1, P1, P1))
fb = np.zeros((P1, P1, P1))
fc = np.zeros((P1, P1, P1))
for n in range(P1):
    for k in range(P1):
        for m in range(0, min(n, k) + 1):
            if k <= n:
                fa[n, k, m] = float(mp.sqrt(comb(n - m, n - k) * comb(n + m, n - k)))
            if k >= n:
                fb[n, k, m] = float(mp.sqrt(comb(k - m, k - n) * comb(k + m, k - n)))
            fc[n, k, m] = float(mp.sqrt(comb(n + k, n + m) * comb(n + k, n - m)))
out += [fa.tobytes(), fb.tobytes(), fc.tobytes()]
nn, mm = H.nm(p)
fe = np.zeros((len(nn), len(nn)))
for i, (n, m) in enumerate(zip(nn, mm)):
    for j, (k, l) in enumerate(zip(nn, mm)):
        fe[i, j] = float(mp.sqrt(comb(n + k - m + l, n - m) * comb(n + k + m - l, n + m)))
out.append(fe.tobytes())
with open(a.o, "wb") as f:
    for b in out:
        f.write(b)
