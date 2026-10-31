# gate4096.py: does |w|max still move beyond 3072? Peaks of Omega/Omega0 and |w|max/|w|0 over t* 1.3-1.95 on the rungs
# l40k1536, l40k2048, l40k2560 and y40k3072 (IC normalization 2.6026367, 20.0). Extrapolation to n -> inf:
#  Aitken (2048, 2560, 3072: equally spaced, exact for spectral convergence), and the fits phi = phi_inf + A n^-p through
#  every three of the four rungs (drop-one-rung); central value = Aitken, spread = max - min over all estimates.
# Moving: |phi_inf - phi(3072)| / phi_inf > 3% for |w|max. Writes yspec/w4096.txt and w4096.go or w4096.no.
import glob, itertools, numpy as np
from scipy.optimize import brentq
Y = "/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec"; TS = 0.4657717102700246; O0, W0 = 2.6026367, 20.0
def rd(fs):
    rows = {}
    for f in fs:
        L = [l.split() for l in open(f) if l.strip()]
        for r in L[1:]:
            rows[int(r[0])] = [float(x) for x in r]
    d = np.array([rows[k] for k in sorted(rows)]); t = d[:, 1] * TS; m = (t >= 1.3) & (t <= 1.95 + 1e-6)
    return t[m], d[m, 3] / O0, d[m, 5] / W0
def peak(t, y):
    i = int(np.argmax(y))
    if 0 < i < len(y) - 1:
        a, b, c = y[i - 1], y[i], y[i + 1]; dd = a - 2 * b + c; x = 0.5 * (a - c) / dd if dd else 0
        return b - 0.25 * (a - c) * x, t[i] + x * (t[i + 1] - t[i])
    return y[i], t[i]
def fit3(ns, fs):
    # phi = phi_inf + A n^-p through three rungs n1 < n2 < n3; None unless monotone and converging
    (n1, n2, n3), (f1, f2, f3) = ns, fs; d1, d2 = f2 - f1, f3 - f2
    if d1 == 0 or d1 * d2 <= 0:
        return None
    g = lambda p: (n2 ** -p - n3 ** -p) / (n1 ** -p - n2 ** -p) - d2 / d1
    try:
        p = brentq(g, 0.05, 40)
    except ValueError:
        return None
    A = d1 / (n2 ** -p - n1 ** -p)
    return f1 - A * n1 ** -p
ns = [1536, 2048, 2560, 3072]
runs = {1536: [Y + "/run/l40k1536/out"], 2048: [Y + "/run/l40k2048/out"], 2560: [Y + "/run/l40k2560/out"],
        3072: sorted(glob.glob(Y + "/run/y40k3072/out.*"))}
P = {n: rd(runs[n]) for n in ns}
L = ["rung   t* end   Omega/Omega0 peak (t*)    |w| ratio peak (t*)"]
pk = {}
for n in ns:
    t, o, w = P[n]; po, to = peak(t, o); pw, tw = peak(t, w); pk[n] = (po, pw)
    L.append("%5d  %.3f    %.4f (%.4f)          %.2f (%.4f)" % (n, t[-1], po, to, pw, tw))
dec = None
for j, name in ((0, "Omega/Omega0 peak"), (1, "|w| ratio peak")):
    f = {n: pk[n][j] for n in ns}
    d1, d2 = f[2560] - f[2048], f[3072] - f[2560]
    ait = f[3072] - d2 * d2 / (d2 - d1) if d2 != d1 and d1 * d2 > 0 and abs(d2) < abs(d1) else None
    est = [e for e in [ait] if e is not None]
    for c in itertools.combinations(ns, 3):
        if 3072 in c:
            e = fit3(c, [f[n] for n in c])
            if e is not None and np.isfinite(e):
                est.append(e)
    if est:
        cen = ait if ait is not None else float(np.median(est)); spr = max(est) - min(est)
        mv = abs(cen - f[3072]) / abs(cen)
        L.append("%s: rungs %s; extrapolated %.4g (Aitken %s), spread %.3g over %d estimates; 3072 is %.2f%% from it"
                 % (name, " / ".join("%.4g" % f[n] for n in ns), cen, "%.4g" % ait if ait is not None else "n/a", spr, len(est), 100 * mv))
    else:
        mv = float("inf")
        L.append("%s: rungs %s; not monotone/converging: no extrapolation" % (name, " / ".join("%.4g" % f[n] for n in ns)))
    if j == 1:
        dec = mv > 0.03
L.append("decision: " + ("GO 4096 window (|w|max moves > 3% beyond 3072)" if dec else "no 4096 window (|w|max within 3%)"))
open(Y + "/w4096.txt", "w").write("\n".join(L) + "\n")
open(Y + ("/w4096.go" if dec else "/w4096.no"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
