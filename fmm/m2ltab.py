import argparse
import sys
import numpy as np
import mpmath as mp

pa = argparse.ArgumentParser(prog="m2ltab")
pa.add_argument("-p", type=int, required=True)
pa.add_argument("-o", required=True)
a = pa.parse_args()
p = a.p
if p < 1 or p > 40:
    sys.stderr.write("m2ltab: error: need 1 <= -p <= 40\n")
    sys.exit(1)
mp.mp.dps = 40
J = 2 * p
nj = (J + 1) * (J + 2) // 2
out = np.zeros((343, nj, 2))
for o in range(343):
    ox, oy, oz = o // 49 - 3, (o // 7) % 7 - 3, o % 7 - 3
    if max(abs(ox), abs(oy), abs(oz)) < 2:
        continue
    x, y, z = mp.mpf(-ox), mp.mpf(-oy), mp.mpf(-oz)
    r2 = x * x + y * y + z * z
    w = mp.mpc(x, y) / r2
    um = 1 / mp.sqrt(r2)
    for m in range(J + 1):
        if m > 0:
            um = um * w * mp.sqrt(mp.mpf(2 * m - 1) / (2 * m))
        u1, u2 = mp.mpc(um), mp.mpc(0)
        out[o, m * (m + 1) // 2 + m] = float(mp.re(u1)), float(mp.im(u1))
        for n in range(m + 1, J + 1):
            u = ((2 * n - 1) * z / r2 * u1 - mp.sqrt((n - 1) ** 2 - m * m) / r2 * u2) / mp.sqrt(n * n - m * m)
            out[o, n * (n + 1) // 2 + m] = float(mp.re(u)), float(mp.im(u))
            u2, u1 = u1, u
with open(a.o, "wb") as f:
    np.array([p], np.int64).tofile(f)
    out.tofile(f)
