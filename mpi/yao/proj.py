# proj.py out C phases...: steps of an adaptive run from the velocity history of a boxsymmp out file (any grid):
# U(t) = (C/dt) / kc(n_ref) = max (|u| + |v| + |w|) (kc = (n - 1) / 3 per direction; out columns t, C, dt; n_ref from
# -n). A phase "n:t0:t1:sps:P" (grid, time window in t*, s/step on P GPUs) costs steps = int kc(n) U / C dt.
# python3 proj.py -n 1024 out.0 1.0 1024:0:1.3:0.5:4 2048:1.3:1.6:1.6:8 3072:1.6:3.0:2.5:16
import sys, numpy as np
TS = 7.680385244 / (2 * np.pi * 1.62 ** 2)
a = sys.argv[1:]
nref = int(a[a.index('-n') + 1]); del a[a.index('-n'):a.index('-n') + 2]
f, C = a[0], float(a[1])
L = [l.split() for l in open(f) if l.strip()]
h = L[0]; d = np.array([[float(x) for x in r] for r in L[1:]]); col = {n: d[:, i] for i, n in enumerate(h)}
t, c, dt = col['t'], col['C'], col['dt']
U = c / dt / ((nref - 1) // 3)
tot = gh = 0
for ph in a[2:]:
    n, t0, t1, sps, P = ph.split(':'); n = int(n); t0, t1 = float(t0) / TS, float(t1) / TS; sps = float(sps); P = int(P)
    tt = np.linspace(t0, t1, 2001)
    if t1 > t[-1]:
        print("  (t* %.3f beyond the reference run: U held at its last value %.3f)" % (t1 * TS, U[-1]))
    Ui = np.interp(tt, t, U)
    steps = np.trapezoid(((n - 1) // 3) * Ui / C, tt)
    w = steps * sps / 3600
    tot += w; gh += w * P
    print("%5d t* %.2f-%.2f: U %.2f-%.2f, dt %.2e-%.2e, steps %7.0f, %.3f s/step on %d GPUs: %6.1f h, %7.0f GPU-h" % (
        n, t0 * TS, t1 * TS, Ui.min(), Ui.max(), C / (((n - 1) // 3) * Ui.max()), C / (((n - 1) // 3) * Ui.min()), steps, sps, P, w, w * P))
print("total %.1f h wall, %.0f GPU-h" % (tot, gh))
