import argparse
import sys
import numpy as np
import mpmath as mp

pa = argparse.ArgumentParser(prog="neartab")
pa.add_argument("-B", type=int, required=True)
pa.add_argument("-K", type=int, default=1)
pa.add_argument("-e", type=float, default=2.0)
pa.add_argument("-o", required=True)
a = pa.parse_args()
if a.B < 1 or a.K < 1:
    sys.stderr.write("neartab: error: need -B >= 1, -K >= 1\n")
    sys.exit(1)
mp.mp.dps = 40
R = (a.K + 1) * a.B - 1
W = 2 * R + 1
dl = mp.mpf(a.e)


def q8(rho):
    x = rho * rho
    return mp.erf(rho) - 2 / mp.sqrt(mp.pi) * rho * mp.exp(-x) * (1 - mp.mpf(89) / 24 * x + mp.mpf(5) / 3 * x * x - x**3 / 6)


full = np.zeros((W, W, W, 3))
diff = np.zeros((W, W, W, 3))
cache = {}
for i in range(-R, R + 1):
    for j in range(-R, R + 1):
        for k in range(-R, R + 1):
            r2 = i * i + j * j + k * k
            if r2 == 0:
                continue
            if r2 not in cache:
                r = mp.sqrt(r2)
                qq = q8(r / dl)
                cache[r2] = (qq / (4 * mp.pi * r**3), (qq - 1) / (4 * mp.pi * r**3))
            f, g = cache[r2]
            for c, v in enumerate((i, j, k)):
                full[i + R, j + R, k + R, c] = float(f * v)
                diff[i + R, j + R, k + R, c] = float(g * v)
with open(a.o, "wb") as fo:
    np.array([a.B, a.K, R], np.int64).tofile(fo)
    np.array([a.e], float).tofile(fo)
    full.tofile(fo)
    diff.tofile(fo)
